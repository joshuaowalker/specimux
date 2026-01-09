"""File watching system for automatic specimux demultiplexing.

This module provides a file watcher that monitors a directory for new FASTQ files
and automatically runs specimux on them. Designed for live MinKNOW sequencing workflows.
"""

import sys
import argparse
import logging
import os
import time
import json
import subprocess
import threading
from pathlib import Path
from typing import Dict, Optional, Set
from datetime import datetime
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler, FileCreatedEvent


class ProcessedFilesTracker:
    """Tracks processed files to avoid reprocessing.

    State is persisted to a JSON file to survive restarts.
    """

    def __init__(self, state_file: Path):
        self.state_file = state_file
        self.processed: Dict[str, dict] = {}
        self._load_state()

    def _load_state(self):
        """Load state from file if it exists."""
        if self.state_file.exists():
            try:
                with open(self.state_file, 'r') as f:
                    data = json.load(f)
                    self.processed = data.get('processed_files', {})
                logging.info(f"Loaded state: {len(self.processed)} previously processed files")
            except Exception as e:
                logging.warning(f"Could not load state file: {e}")
                self.processed = {}

    def _save_state(self):
        """Save current state to file."""
        try:
            with open(self.state_file, 'w') as f:
                json.dump({'processed_files': self.processed}, f, indent=2)
        except Exception as e:
            logging.error(f"Could not save state file: {e}")

    def is_processed(self, filepath: str) -> bool:
        """Check if file has been processed."""
        return filepath in self.processed

    def mark_processed(self, filepath: str, status: str, size: int):
        """Mark file as processed with status (success/failed)."""
        self.processed[filepath] = {
            'timestamp': datetime.now().isoformat(),
            'size': size,
            'status': status
        }
        self._save_state()

    def get_status(self, filepath: str) -> Optional[str]:
        """Get processing status for a file."""
        return self.processed.get(filepath, {}).get('status')


class FileStabilityChecker:
    """Checks if a file has finished being written by monitoring size stability."""

    def __init__(self, settle_time: int = 30):
        self.settle_time = settle_time

    def wait_for_stable(self, filepath: Path, check_interval: int = 5) -> bool:
        """Wait until file size is stable for settle_time seconds.

        Returns True if stable, False if file disappeared or check_interval exceeded settle_time.

        Args:
            filepath: Path to file to monitor
            check_interval: Seconds between size checks (default: 5)
        """
        if not filepath.exists():
            logging.warning(f"File {filepath} does not exist")
            return False

        try:
            last_size = filepath.stat().st_size
            stable_duration = 0

            logging.info(f"Waiting for {filepath.name} to stabilize (initial size: {last_size} bytes)")

            while stable_duration < self.settle_time:
                time.sleep(check_interval)

                if not filepath.exists():
                    logging.warning(f"File {filepath} disappeared during stability check")
                    return False

                current_size = filepath.stat().st_size

                if current_size != last_size:
                    # Size changed, reset timer
                    logging.debug(f"{filepath.name} size changed: {last_size} -> {current_size} bytes")
                    last_size = current_size
                    stable_duration = 0
                else:
                    # Size unchanged, increment stable duration
                    stable_duration += check_interval
                    logging.debug(f"{filepath.name} stable for {stable_duration}s (target: {self.settle_time}s)")

            logging.info(f"{filepath.name} is stable at {last_size} bytes")
            return True

        except Exception as e:
            logging.error(f"Error checking file stability: {e}")
            return False


class SpecimuxRunner:
    """Executes specimux on stable files."""

    def __init__(self, primers_file: str, specimens_file: str, specimux_args: list):
        self.primers_file = primers_file
        self.specimens_file = specimens_file
        self.specimux_args = specimux_args

    def run(self, fastq_file: Path) -> tuple[bool, str]:
        """Run specimux on the given file.

        Returns:
            (success: bool, stderr_output: str)
        """
        cmd = [
            'specimux',
            self.primers_file,
            self.specimens_file,
            str(fastq_file)
        ] + self.specimux_args

        logging.info(f"Running: {' '.join(cmd)}")

        try:
            # Don't capture stdout or stderr - let everything pass through
            # This allows progress bars and logging to display in real-time
            result = subprocess.run(
                cmd,
                stdout=None,  # Inherit stdout
                stderr=None,  # Inherit stderr (progress bars use this)
                timeout=None  # No timeout for potentially long runs
            )

            if result.returncode == 0:
                logging.info(f"Successfully processed {fastq_file.name}")
                return True, ""
            else:
                logging.error(f"Specimux failed on {fastq_file.name} (exit code: {result.returncode})")
                return False, f"Exit code: {result.returncode}"

        except subprocess.TimeoutExpired:
            logging.error(f"Specimux timed out on {fastq_file.name}")
            return False, "Timeout"
        except Exception as e:
            logging.error(f"Error running specimux on {fastq_file.name}: {e}")
            return False, str(e)


class FastqFileHandler(FileSystemEventHandler):
    """Handles file system events for new FASTQ files.

    Uses a lock to ensure only one file is processed at a time, even if
    multiple files are detected simultaneously.
    """

    def __init__(
        self,
        runner: SpecimuxRunner,
        tracker: ProcessedFilesTracker,
        stability_checker: FileStabilityChecker,
        pattern: str = "*.fastq",
        stop_after: Optional[int] = None
    ):
        self.runner = runner
        self.tracker = tracker
        self.stability_checker = stability_checker
        self.pattern = pattern
        self.stop_after = stop_after
        self.processed_count = 0
        self.pending_files: Set[str] = set()
        self.processing_lock = threading.Lock()  # Ensure sequential processing

    def on_created(self, event):
        """Handle file creation events.

        Note: This method may be called from multiple threads, but the processing_lock
        ensures only one file is processed at a time.
        """
        if event.is_directory:
            return

        filepath = Path(event.src_path)

        # Check if file matches pattern
        if not filepath.match(self.pattern):
            return

        filename = str(filepath)

        # Skip if already processed (quick check without lock)
        if self.tracker.is_processed(filename):
            status = self.tracker.get_status(filename)
            logging.info(f"Skipping {filepath.name} (already processed: {status})")
            return

        # Acquire lock to ensure sequential processing
        # If another file is being processed, this will block until it's done
        if self.processing_lock.locked():
            logging.info(f"File {filepath.name} detected, waiting for current processing to finish...")

        with self.processing_lock:
            # Check again after acquiring lock (another thread might have started)
            if filename in self.pending_files:
                logging.debug(f"Skipping {filepath.name} (already pending)")
                return

            if self.tracker.is_processed(filename):
                status = self.tracker.get_status(filename)
                logging.debug(f"Skipping {filepath.name} (processed while waiting)")
                return

            logging.info(f"New file detected: {filepath.name}")
            self.pending_files.add(filename)

            # Process the file (still holding lock)
            self._process_file(filepath)

    def _process_file(self, filepath: Path):
        """Process a single FASTQ file."""
        filename = str(filepath)

        try:
            # Wait for file to be stable
            if not self.stability_checker.wait_for_stable(filepath):
                logging.warning(f"File {filepath.name} did not stabilize, skipping")
                self.tracker.mark_processed(filename, 'failed', 0)
                return

            # Get final file size
            file_size = filepath.stat().st_size

            # Run specimux
            success, output = self.runner.run(filepath)

            # Mark as processed
            status = 'success' if success else 'failed'
            self.tracker.mark_processed(filename, status, file_size)

            self.processed_count += 1
            logging.info(f"Processed {self.processed_count} file(s) total")

            # Check if we should stop
            if self.stop_after and self.processed_count >= self.stop_after:
                logging.info(f"Reached stop limit of {self.stop_after} files")
                # Signal to stop (observer will check this)

        except Exception as e:
            logging.error(f"Error processing {filepath.name}: {e}")
            self.tracker.mark_processed(filename, 'failed', 0)
        finally:
            self.pending_files.discard(filename)


def setup_logging(daemon: bool, output_dir: Optional[str] = None, debug: bool = False):
    """Set up logging for the watcher.

    In normal mode: log to stdout
    In daemon mode: log to file in output_dir or current directory
    """
    level = logging.DEBUG if debug else logging.INFO
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')

    logging.getLogger().handlers.clear()

    if daemon:
        # Daemon mode: log to file
        log_dir = Path(output_dir) if output_dir else Path.cwd()
        log_dir.mkdir(parents=True, exist_ok=True)
        log_file = log_dir / 'specimux-watch.log'

        file_handler = logging.FileHandler(log_file, mode='a')
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)

        print(f"Daemon mode: logging to {log_file}")
    else:
        # Normal mode: log to console
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logging.getLogger().addHandler(console_handler)

    logging.getLogger().setLevel(level)


def parse_args(argv):
    """Parse command-line arguments for specimux-watch."""
    parser = argparse.ArgumentParser(
        description="Watch directory for new FASTQ files and automatically run specimux.",
        epilog="Example: specimux-watch primers.fasta specimens.txt /path/to/watch -F -O output/ -d"
    )

    # Required arguments (matching specimux order)
    parser.add_argument("primer_file", help="Fasta file containing primer information")
    parser.add_argument("specimen_file", help="TSV file containing specimen mapping with barcodes and primers")
    parser.add_argument("watch_dir", help="Directory to watch for new FASTQ files")

    # Watch-specific options
    parser.add_argument("--settle-time", type=int, default=30,
                       help="Seconds to wait for file size stability (default: 30)")
    parser.add_argument("--state-file", type=str, default=None,
                       help="Path to state file for tracking processed files (default: .specimux-watch-state.json in watch dir)")
    parser.add_argument("--pattern", type=str, default="*.fastq",
                       help="File pattern to watch (default: *.fastq)")
    parser.add_argument("--daemon", action="store_true",
                       help="Run in daemon mode (log to file instead of stdout)")
    parser.add_argument("--stop-after", type=int, default=None,
                       help="Stop after processing N files (useful for testing)")

    # Common specimux arguments that should be passed through
    parser.add_argument("--min-length", type=int, default=-1)
    parser.add_argument("--max-length", type=int, default=-1)
    parser.add_argument("-e", "--index-edit-distance", type=int, default=-1)
    parser.add_argument("-E", "--primer-edit-distance", type=int, default=-1)
    parser.add_argument("-l", "--search-len", type=int, default=80)
    parser.add_argument("-F", "--output-to-files", action="store_true")
    parser.add_argument("-P", "--output-file-prefix", default="")
    parser.add_argument("-O", "--output-dir", default=".")
    parser.add_argument("--color", action="store_true")
    parser.add_argument("--trim", type=str, default="barcodes")
    parser.add_argument("--dereplicate", type=str, default="best")
    parser.add_argument("-d", "--diagnostics", nargs='?', const=1, type=int, choices=[1, 2, 3])
    parser.add_argument("-D", "--debug", action="store_true")
    parser.add_argument("--disable-prefilter", action="store_true")
    parser.add_argument("--disable-preorient", action="store_true")
    parser.add_argument("-t", "--threads", type=int, default=-1)
    parser.add_argument("--sample-topq", type=int, default=0)

    args = parser.parse_args(argv[1:])

    # Validate watch directory exists
    watch_path = Path(args.watch_dir)
    if not watch_path.exists():
        parser.error(f"Watch directory does not exist: {args.watch_dir}")
    if not watch_path.is_dir():
        parser.error(f"Watch path is not a directory: {args.watch_dir}")

    # Set default state file if not specified
    if args.state_file is None:
        args.state_file = str(watch_path / '.specimux-watch-state.json')

    return args


def build_specimux_args(args) -> list:
    """Build the list of arguments to pass to specimux."""
    specimux_args = []

    if args.min_length != -1:
        specimux_args.extend(['--min-length', str(args.min_length)])
    if args.max_length != -1:
        specimux_args.extend(['--max-length', str(args.max_length)])
    if args.index_edit_distance != -1:
        specimux_args.extend(['-e', str(args.index_edit_distance)])
    if args.primer_edit_distance != -1:
        specimux_args.extend(['-E', str(args.primer_edit_distance)])
    if args.search_len != 80:
        specimux_args.extend(['-l', str(args.search_len)])
    if args.output_to_files:
        specimux_args.append('-F')
    if args.output_file_prefix != "":
        specimux_args.extend(['-P', args.output_file_prefix])
    if args.output_dir != ".":
        specimux_args.extend(['-O', args.output_dir])
    if args.color:
        specimux_args.append('--color')
    if args.trim != "barcodes":
        specimux_args.extend(['--trim', args.trim])
    if args.dereplicate != "best":
        specimux_args.extend(['--dereplicate', args.dereplicate])
    if args.diagnostics:
        specimux_args.extend(['-d', str(args.diagnostics)])
    if args.debug:
        specimux_args.append('-D')
    if args.disable_prefilter:
        specimux_args.append('--disable-prefilter')
    if args.disable_preorient:
        specimux_args.append('--disable-preorient')
    if args.threads != -1:
        specimux_args.extend(['-t', str(args.threads)])
    if args.sample_topq != 0:
        specimux_args.extend(['--sample-topq', str(args.sample_topq)])

    return specimux_args


def main():
    """Main entry point for specimux-watch command."""
    args = parse_args(sys.argv)

    # Set up logging
    setup_logging(args.daemon, args.output_dir if args.output_to_files else None, args.debug)

    logging.info("Starting specimux-watch")
    logging.info(f"Watching directory: {args.watch_dir}")
    logging.info(f"Pattern: {args.pattern}")
    logging.info(f"Settle time: {args.settle_time}s")
    logging.info(f"State file: {args.state_file}")

    # Build specimux arguments
    specimux_args = build_specimux_args(args)
    logging.info(f"Specimux arguments: {' '.join(specimux_args)}")

    # Initialize components with fresh state
    state_file = Path(args.state_file)

    # Remove old state file if it exists
    if state_file.exists():
        logging.info("Removing old state file")
        state_file.unlink()

    tracker = ProcessedFilesTracker(state_file)

    # Mark all pre-existing files as already processed (ignored)
    watch_path = Path(args.watch_dir)
    existing_files = list(watch_path.glob(args.pattern))
    if existing_files:
        logging.info(f"Found {len(existing_files)} pre-existing file(s) - marking as ignored")
        for filepath in existing_files:
            try:
                size = filepath.stat().st_size
                tracker.mark_processed(str(filepath), 'ignored', size)
                logging.debug(f"  Ignoring: {filepath.name}")
            except Exception as e:
                logging.warning(f"Could not stat {filepath.name}: {e}")

    stability_checker = FileStabilityChecker(settle_time=args.settle_time)
    runner = SpecimuxRunner(args.primer_file, args.specimen_file, specimux_args)

    # Set up file handler and observer
    event_handler = FastqFileHandler(
        runner=runner,
        tracker=tracker,
        stability_checker=stability_checker,
        pattern=args.pattern,
        stop_after=args.stop_after
    )

    observer = Observer()
    observer.schedule(event_handler, args.watch_dir, recursive=False)
    observer.start()

    logging.info("Watching for new files (Ctrl+C to stop)...")

    try:
        while True:
            time.sleep(1)
            # Check if we should stop
            if args.stop_after and event_handler.processed_count >= args.stop_after:
                logging.info("Stopping after processing requested number of files")
                break
    except KeyboardInterrupt:
        logging.info("Received interrupt, stopping...")
    finally:
        observer.stop()
        observer.join()
        logging.info("Stopped")


if __name__ == "__main__":
    main()
