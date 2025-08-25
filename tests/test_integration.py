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
        
        assert hasattr(MultipleMatchStrategy, 'RETAIN')
        assert hasattr(MultipleMatchStrategy, 'DOWNGRADE_FULL')


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