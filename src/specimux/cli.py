"""Command-line interface for specimux tools."""
import sys
import argparse
import logging
import os
from . import core, __version__
from .core import TrimMode, MultipleMatchStrategy


def version():
    """Return version string."""
    return f"specimux version {__version__}"


def parse_args(argv):
    """Parse command-line arguments for specimux."""
    parser = argparse.ArgumentParser(description="Specimux: Demultiplex MinION sequences by dual barcode indexes and primers.")

    parser.add_argument("primer_file", help="Fasta file containing primer information")
    parser.add_argument("specimen_file", help="TSV file containing specimen mapping with barcodes and primers")
    parser.add_argument("sequence_file", help="Sequence file in Fasta or Fastq format, gzipped or plain text")

    parser.add_argument("--min-length", type=int, default=-1, help="Minimum sequence length.  Shorter sequences will be skipped (default: no filtering)")
    parser.add_argument("--max-length", type=int, default=-1, help="Maximum sequence length.  Longer sequences will be skipped (default: no filtering)")
    parser.add_argument("-n", "--num-seqs", type=str, default="-1", help="Number of sequences to read from file (e.g., -n 100 or -n 102,3)")
    parser.add_argument("-e", "--index-edit-distance", type=int, default=-1, help="Barcode edit distance value, default is half of min distance between barcodes")
    parser.add_argument("-E", "--primer-edit-distance", type=int, default=-1, help="Primer edit distance value, default is min distance between primers")
    parser.add_argument("-l", "--search-len", type=int, default=80, help="Length to search for index and primer at start and end of sequence (default: 80)")
    parser.add_argument("-F", "--output-to-files", action="store_true", help="Create individual sample files for sequences")
    parser.add_argument("-P", "--output-file-prefix", default="", help="Prefix for individual files when using -F (default: no prefix)")
    parser.add_argument("-O", "--output-dir", default=".", help="Directory for individual files when using -F (default: .)")
    parser.add_argument("--color", action="store_true", help="Highlight barcode matches in blue, primer matches in green")
    parser.add_argument("--trim", choices=[TrimMode.NONE, TrimMode.TAILS, TrimMode.BARCODES, TrimMode.PRIMERS], default=TrimMode.BARCODES, help="trimming to apply")
    parser.add_argument("--dereplicate", choices=[MultipleMatchStrategy.NONE, MultipleMatchStrategy.BEST], default=MultipleMatchStrategy.BEST, help="Dereplication strategy: 'best' selects best match per specimen/barcode group (default), 'none' outputs all matches")
    parser.add_argument("-d", "--diagnostics", nargs='?', const=1, type=int, choices=[1, 2, 3], 
                        help="Enable diagnostic trace logging: 1=standard (default), 2=detailed, 3=verbose")
    parser.add_argument("-D", "--debug", action="store_true", help="Enable debug logging")

    parser.add_argument("--disable-prefilter", action="store_true", help="Disable barcode prefiltering (bloom filter optimization)")
    parser.add_argument("--disable-preorient", action="store_true", help="Disable heuristic pre-orientation")
    parser.add_argument("-t", "--threads", type=int, default=-1, help="Number of worker threads to use")
    parser.add_argument("--sample-topq", type=int, default=0, metavar="N",
                        help="Create subsample directories with top N sequences by average quality score (default: disabled)")
    parser.add_argument("-v", "--version", action="version", version=version())

    args = parser.parse_args(argv[1:])

    if args.num_seqs:
        process_num_seqs(args, parser)

    return args


def process_num_seqs(args, parser):
    """Process the -n/--num-seqs argument."""
    if ',' in args.num_seqs:
        start, num = args.num_seqs.split(',')
        try:
            args.start_seq = int(start)
            args.num_seqs = int(num)
        except ValueError:
            parser.error("Invalid format for -n option. Use 'start,num' with integers.")
    else:
        try:
            args.num_seqs = int(args.num_seqs)
            args.start_seq = 1
        except ValueError:
            parser.error("Invalid format for -n option. Use an integer or 'start,num' with integers.")


def setup_logging(debug: bool, output_dir: str = None, is_worker: bool = False):
    """Set up logging to both console and file if output directory is specified."""
    level = logging.DEBUG if debug else logging.INFO
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Clear any existing handlers
    logging.getLogger().handlers.clear()
    
    # Set up console logging
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logging.getLogger().addHandler(console_handler)
    
    # Set up file logging if output directory specified
    # Only the main process should create the file handler
    if output_dir and not is_worker:
        os.makedirs(output_dir, exist_ok=True)  # Ensure directory exists
        log_file = os.path.join(output_dir, 'log.txt')
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setFormatter(formatter)
        logging.getLogger().addHandler(file_handler)
        
    logging.getLogger().setLevel(level)


def main():
    """Main entry point for specimux command."""
    args = parse_args(sys.argv)
    setup_logging(args.debug, args.output_dir if args.output_to_files else None)
    
    # Log version and command line used
    logging.info(f"Starting {version()}")
    logging.info(f"Command line: {' '.join(sys.argv)}")

    if args.output_to_files:
        core.specimux_mp(args)  # Use multiprocess for file output
    else:
        if args.threads > 1:
            logging.warning(f"Multithreading only supported for file output. Ignoring --threads {args.threads}")
        core.specimux(args)     # Use single process for console output


def specimine_main():
    """Entry point for specimine command."""
    from . import specimine
    specimine.main()


def converter_main():
    """Entry point for specimux-convert command."""
    from . import converter
    converter.main()


def trace_main():
    """Entry point for specimux-stats command."""
    from . import trace_stats
    trace_stats.main()


def visualize_main():
    """Entry point for specimux-visualize command."""
    from . import visualize
    visualize.main()


def watch_main():
    """Entry point for specimux-watch command."""
    from . import watch
    watch.main()


if __name__ == "__main__":
    main()