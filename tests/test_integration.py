#!/usr/bin/env python3
"""
Integration tests for specimux using pytest.
Tests the complete pipeline with real data.
"""

import subprocess
import sys
import tempfile
import shutil
from pathlib import Path
import pytest
import json


class TestSpecimuxIntegration:
    """Integration tests for the specimux pipeline."""
    
    @staticmethod
    def extract_match_rate(stderr_output):
        """Extract match rate from specimux stderr output.
        
        Args:
            stderr_output: The stderr text from specimux run
            
        Returns:
            float: The match rate percentage
            
        Raises:
            AssertionError: If match rate cannot be found
        """
        import re
        match = re.search(r"match rate: ([\d.]+)%", stderr_output)
        assert match, f"Could not find match rate in output: {stderr_output}"
        return float(match.group(1))
    
    @staticmethod
    def assert_match_rate_in_range(actual_rate, expected_rate, tolerance=5.0):
        """Assert that match rate is within acceptable range.
        
        Args:
            actual_rate: The actual match rate from the test
            expected_rate: The expected baseline match rate
            tolerance: Acceptable variance in percentage points (default 5%)
                      This accounts for small sample size variations and
                      ensures the algorithm hasn't degraded.
        """
        assert abs(actual_rate - expected_rate) <= tolerance, \
            f"Match rate {actual_rate}% outside expected range " \
            f"[{expected_rate - tolerance}, {expected_rate + tolerance}]. " \
            f"This may indicate algorithm degradation."
    
    @pytest.fixture(scope="class")
    def test_data_dir(self):
        """Fixture providing path to test data directory."""
        return Path(__file__).parent / "data" / "integration_test_suite"
    
    @pytest.fixture(scope="class")
    def test_files(self, test_data_dir):
        """Fixture providing paths to test input files."""
        return {
            'primers': test_data_dir / "primers.fasta",
            'specimens': test_data_dir / "specimens.txt",
            'sequences': test_data_dir / "sequences.fastq",
            'expected_output': test_data_dir / "expected_output",
            'validation_script': test_data_dir / "validate_test_results.py"
        }
    
    @pytest.fixture
    def temp_output_dir(self):
        """Fixture providing a temporary output directory."""
        with tempfile.TemporaryDirectory(prefix="specimux_test_") as temp_dir:
            yield Path(temp_dir) / "output"
    
    @pytest.mark.integration
    def test_full_pipeline(self, test_files, temp_output_dir):
        """Test the complete specimux pipeline with 40 test sequences."""
        # Verify test files exist
        for name, path in test_files.items():
            assert path.exists(), f"Required test file missing: {name} at {path}"
        
        # Run specimux
        cmd = [
            sys.executable, "-m", "specimux.cli",
            str(test_files['primers']),
            str(test_files['specimens']),
            str(test_files['sequences']),
            "-F", "-O", str(temp_output_dir), "-d"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Specimux failed: {result.stderr}"
        # Check in stderr since logging goes there
        assert "Processed 40 sequences" in result.stderr
        
        # Verify match rate is within expected range (15% ± 5% for this test dataset)
        actual_rate = self.extract_match_rate(result.stderr)
        expected_rate = 15.0  # Historical baseline for this test dataset
        self.assert_match_rate_in_range(actual_rate, expected_rate)
        
        # Validate output structure
        self._validate_output_structure(temp_output_dir, test_files['expected_output'])
        
        # Run validation script
        validation_cmd = [
            sys.executable,
            str(test_files['validation_script']),
            str(temp_output_dir),
            str(test_files['expected_output'])
        ]
        
        validation_result = subprocess.run(validation_cmd, capture_output=True, text=True)
        assert validation_result.returncode == 0, f"Validation failed: {validation_result.stdout}"
        assert "All tests PASSED" in validation_result.stdout
    
    @pytest.mark.integration
    @pytest.mark.parametrize("num_sequences,expected_match_rate", [
        (5, 20.0),   # First 5 sequences - 1 match
        (10, 20.0),  # First 10 sequences - 2 matches
        (20, 20.0),  # First 20 sequences - 4 matches
    ])
    def test_partial_sequences(self, test_files, temp_output_dir, num_sequences, expected_match_rate):
        """Test specimux with different numbers of sequences."""
        cmd = [
            sys.executable, "-m", "specimux.cli",
            str(test_files['primers']),
            str(test_files['specimens']),
            str(test_files['sequences']),
            "-n", str(num_sequences),
            "-F", "-O", str(temp_output_dir)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Specimux failed: {result.stderr}"
        # Check in stderr since logging goes there
        assert f"Processed {num_sequences} sequences" in result.stderr
        
        # Extract and verify match rate within acceptable range
        actual_rate = self.extract_match_rate(result.stderr)
        self.assert_match_rate_in_range(actual_rate, expected_match_rate)
    
    def _validate_output_structure(self, actual_dir, expected_dir):
        """Validate that output directory structure matches expected."""
        # Check key directories exist
        assert (actual_dir / "full").exists(), "Missing 'full' directory"
        assert (actual_dir / "partial").exists(), "Missing 'partial' directory"
        assert (actual_dir / "unknown").exists(), "Missing 'unknown' directory"
        assert (actual_dir / "trace").exists(), "Missing 'trace' directory"
        
        # Check for expected files in full matches
        full_its2 = actual_dir / "full" / "ITS2"
        assert full_its2.exists(), "Missing full/ITS2 directory"
        
        # Count files in key locations
        full_files = list((actual_dir / "full").rglob("*.fastq"))
        assert len(full_files) >= 2, f"Expected at least 2 full match files, got {len(full_files)}"


class TestSpecimuxCommands:
    """Test individual specimux commands."""
    
    @pytest.mark.unit
    def test_specimux_version(self):
        """Test specimux version command."""
        result = subprocess.run(
            [sys.executable, "-m", "specimux.cli", "--version"],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        assert "specimux version" in result.stdout
    
    @pytest.mark.unit
    def test_specimux_help(self):
        """Test specimux help command."""
        result = subprocess.run(
            [sys.executable, "-m", "specimux.cli", "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        assert "Demultiplex MinION sequences" in result.stdout
        assert "primer_file" in result.stdout
        assert "specimen_file" in result.stdout
    
    @pytest.mark.unit
    def test_specimine_help(self):
        """Test specimine help command."""
        from specimux.cli import specimine_main
        
        # Test that the function exists and is callable
        assert callable(specimine_main)
        
        # Test command line help
        result = subprocess.run(
            [sys.executable, "-c", 
             "from specimux.cli import specimine_main; import sys; sys.argv=['specimine', '--help']; specimine_main()"],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        assert "Mine additional candidate sequences" in result.stdout


class TestCoreModules:
    """Unit tests for core specimux modules."""
    
    @pytest.mark.unit
    def test_primer_database_import(self):
        """Test that PrimerDatabase can be imported and instantiated."""
        from specimux import PrimerDatabase
        
        # Create empty primer database
        db = PrimerDatabase()
        assert db is not None
        # Check for expected methods
        assert hasattr(db, 'add_primer')
        assert hasattr(db, 'get_primer')
    
    @pytest.mark.unit
    def test_specimens_import(self):
        """Test that Specimens class can be imported."""
        from specimux import Specimens
        
        assert Specimens is not None
        assert hasattr(Specimens, '__init__')
    
    @pytest.mark.unit
    def test_match_parameters_import(self):
        """Test that MatchParameters can be imported."""
        from specimux import MatchParameters
        
        assert MatchParameters is not None
    
    @pytest.mark.unit
    def test_trim_modes(self):
        """Test that TrimMode enum is accessible."""
        from specimux.core import TrimMode
        
        assert hasattr(TrimMode, 'NONE')
        assert hasattr(TrimMode, 'TAILS')
        assert hasattr(TrimMode, 'BARCODES')
        assert hasattr(TrimMode, 'PRIMERS')
    
    @pytest.mark.unit
    def test_multiple_match_strategy(self):
        """Test that MultipleMatchStrategy enum is accessible."""
        from specimux.core import MultipleMatchStrategy

        assert hasattr(MultipleMatchStrategy, 'NONE')
        assert hasattr(MultipleMatchStrategy, 'BEST')


class TestSpecimuxWatch:
    """Integration tests for specimux-watch file monitoring."""

    @pytest.fixture
    def watch_setup(self):
        """Fixture providing temporary directories and test files for watch tests."""
        with tempfile.TemporaryDirectory(prefix="specimux_watch_test_") as temp_dir:
            temp_path = Path(temp_dir)
            watch_dir = temp_path / "watch"
            output_dir = temp_path / "output"
            watch_dir.mkdir()

            # Get test data paths
            test_data_dir = Path(__file__).parent / "data" / "integration_test_suite"

            yield {
                'watch_dir': watch_dir,
                'output_dir': output_dir,
                'test_fastq': test_data_dir / "sequences.fastq",
                'primers': test_data_dir / "primers.fasta",
                'specimens': test_data_dir / "specimens.txt"
            }

    @pytest.mark.integration
    def test_watch_ignores_existing_and_processes_new(self, watch_setup):
        """Test that specimux-watch ignores pre-existing files and processes new ones."""
        import time

        # Copy a pre-existing file
        existing_file = watch_setup['watch_dir'] / "existing.fastq"
        shutil.copy(watch_setup['test_fastq'], existing_file)

        # Start specimux-watch in background
        cmd = [
            sys.executable, "-m", "specimux.watch",
            str(watch_setup['primers']),
            str(watch_setup['specimens']),
            str(watch_setup['watch_dir']),
            "-F", "-O", str(watch_setup['output_dir']),
            "--settle-time", "2",
            "--stop-after", "1"  # Stop after processing 1 file
        ]

        watch_process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Give watcher time to start and scan for existing files
        time.sleep(1)

        # Copy a new file
        new_file = watch_setup['watch_dir'] / "newfile.fastq"
        shutil.copy(watch_setup['test_fastq'], new_file)

        # Wait for watch process to complete
        stdout, stderr = watch_process.communicate(timeout=30)
        combined_output = stdout + stderr

        # Check that it ignored existing file
        assert "pre-existing file(s)" in combined_output
        assert "marking as ignored" in combined_output

        # Check that it processed new file
        assert "New file detected: newfile.fastq" in combined_output
        assert "Successfully processed newfile.fastq" in combined_output or \
               "Processed 1 file(s) total" in combined_output

        # Verify state file exists and has correct structure
        state_file = watch_setup['watch_dir'] / ".specimux-watch-state.json"
        assert state_file.exists(), "State file was not created"

        with open(state_file) as f:
            state = json.load(f)

        # Check that existing file was marked as ignored
        existing_entries = [
            entry for path, entry in state['processed_files'].items()
            if 'existing.fastq' in path
        ]
        assert len(existing_entries) == 1
        assert existing_entries[0]['status'] == 'ignored'

        # Check that new file was processed
        new_entries = [
            entry for path, entry in state['processed_files'].items()
            if 'newfile.fastq' in path
        ]
        assert len(new_entries) == 1
        assert new_entries[0]['status'] in ['success', 'failed']  # Either is ok for this test

    @pytest.mark.integration
    def test_watch_help(self):
        """Test that specimux-watch help command works."""
        result = subprocess.run(
            [sys.executable, "-m", "specimux.watch", "--help"],
            capture_output=True,
            text=True
        )
        assert result.returncode == 0
        assert "Watch directory for new FASTQ files" in result.stdout
        assert "watch_dir" in result.stdout
        assert "--settle-time" in result.stdout


@pytest.mark.slow
@pytest.mark.integration
class TestLargeDataset:
    """Tests with larger datasets (marked as slow)."""

    @pytest.fixture(scope="class")
    def large_dataset_path(self):
        """Path to large test dataset (if available)."""
        # This would point to ont37 or similar large dataset
        path = Path("/path/to/large/dataset")
        if not path.exists():
            pytest.skip("Large dataset not available")
        return path

    def test_performance_with_large_dataset(self, large_dataset_path):
        """Test performance with large dataset."""
        # This is a placeholder for performance testing
        # Would run specimux on larger dataset and check timing
        pass


if __name__ == "__main__":
    # Allow running as a script for backwards compatibility
    pytest.main([__file__, "-v"])