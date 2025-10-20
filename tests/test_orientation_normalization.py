#!/usr/bin/env python3
"""
Integration test to verify sequence orientation normalization in specimux.

This test verifies that sequences are output in a normalized orientation
regardless of their input orientation. Both forward and reverse complement
sequences of the same biological sample should produce identical output sequences.
"""

import pytest
import tempfile
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import os


class TestOrientationNormalization:
    """Integration tests for sequence orientation normalization behavior."""

    @pytest.fixture
    def test_data_dir(self):
        """Get path to test data directory."""
        return Path(__file__).parent / "data" / "integration_test_suite"

    def test_orientation_normalization(self, test_data_dir):
        """
        Test that sequences are normalized to the same orientation in output.
        
        This integration test:
        1. Runs specimux on original test sequences
        2. Runs specimux on reverse complement test sequences  
        3. Verifies that matching sequences are output in the same normalized orientation
           (i.e., both forward and reverse input sequences produce identical output sequences)
        """
        
        # Create temporary directories for outputs
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            original_output = temp_path / "original_output"
            rc_output = temp_path / "rc_output"
            
            # Files to use
            primers_file = test_data_dir / "primers.fasta"
            specimens_file = test_data_dir / "specimens.txt"
            original_sequences = test_data_dir / "sequences.fastq"
            rc_sequences = test_data_dir / "sequences_rc.fastq"
            
            # Run specimux on original sequences
            cmd_original = [
                "python", "-m", "specimux.cli",
                str(primers_file),
                str(specimens_file), 
                str(original_sequences),
                "-F", "-O", str(original_output)
            ]
            
            result_original = subprocess.run(
                cmd_original,
                capture_output=True,
                text=True,
                cwd=Path(__file__).parent.parent  # Run from specimux root
            )
            
            assert result_original.returncode == 0, f"Original specimux run failed: {result_original.stderr}"
            
            # Run specimux on reverse complement sequences
            cmd_rc = [
                "python", "-m", "specimux.cli",
                str(primers_file),
                str(specimens_file),
                str(rc_sequences),
                "-F", "-O", str(rc_output)
            ]
            
            result_rc = subprocess.run(
                cmd_rc,
                capture_output=True,
                text=True,
                cwd=Path(__file__).parent.parent  # Run from specimux root
            )
            
            assert result_rc.returncode == 0, f"RC specimux run failed: {result_rc.stderr}"
            
            # Find successful matches in both outputs
            original_full_dir = original_output / "full"
            rc_full_dir = rc_output / "full"
            
            if not original_full_dir.exists():
                pytest.skip("No successful matches found in original sequences - cannot test normalization")
            
            if not rc_full_dir.exists():
                pytest.skip("No successful matches found in RC sequences - cannot test normalization")
            
            # Compare corresponding output files - they should be IDENTICAL
            # because orientation should be normalized
            normalization_matches_found = 0
            
            for pool_dir in original_full_dir.iterdir():
                if pool_dir.is_dir():
                    rc_pool_dir = rc_full_dir / pool_dir.name
                    if rc_pool_dir.exists():
                        # Check pool-level files
                        for original_file in pool_dir.glob("*.fastq"):
                            rc_file = rc_pool_dir / original_file.name
                            if rc_file.exists():
                                # Compare sequences in these files
                                original_seqs = list(SeqIO.parse(original_file, "fastq"))
                                rc_seqs = list(SeqIO.parse(rc_file, "fastq"))
                                
                                if len(original_seqs) > 0 and len(rc_seqs) > 0:
                                    # Take first sequence from each as a representative
                                    orig_seq = str(original_seqs[0].seq)
                                    rc_seq = str(rc_seqs[0].seq)
                                    
                                    # CORRECT BEHAVIOR: These should be IDENTICAL
                                    # because orientation should be normalized
                                    assert orig_seq == rc_seq, (
                                        f"Orientation normalization failed in {original_file.name}!\n"
                                        f"Original output: {orig_seq[:50]}...\n"
                                        f"RC output:       {rc_seq[:50]}...\n"
                                        f"Expected: Identical sequences (normalized orientation)\n"
                                        f"Actual: Different sequences - orientation not normalized\n"
                                        f"This indicates specimux is not properly normalizing sequence orientation."
                                    )
                                    
                                    normalization_matches_found += 1
                                    
                                    # Only need to verify a few examples
                                    if normalization_matches_found >= 3:
                                        break
            
            # Ensure we actually found some matches to test
            assert normalization_matches_found > 0, (
                "No matching specimens found between original and RC runs to test normalization behavior. "
                "This could indicate a problem with the test setup or the sequences don't match the specimens."
            )
            
            print(f"Successfully verified orientation normalization for {normalization_matches_found} specimen matches")

    def test_both_orientations_match_same_specimens(self, test_data_dir):
        """
        Test that forward and reverse sequences of the same biological sample 
        get assigned to the same specimen ID.
        
        This ensures the matching logic works correctly in both orientations.
        """
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            original_output = temp_path / "original_output" 
            rc_output = temp_path / "rc_output"
            
            # Files to use
            primers_file = test_data_dir / "primers.fasta"
            specimens_file = test_data_dir / "specimens.txt"
            original_sequences = test_data_dir / "sequences.fastq"
            rc_sequences = test_data_dir / "sequences_rc.fastq"
            
            # Run specimux on both sequence sets
            for sequences, output_dir in [(original_sequences, original_output), (rc_sequences, rc_output)]:
                cmd = [
                    "python", "-m", "specimux.cli",
                    str(primers_file),
                    str(specimens_file),
                    str(sequences), 
                    "-F", "-O", str(output_dir)
                ]
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    cwd=Path(__file__).parent.parent
                )
                
                assert result.returncode == 0, f"Specimux run failed: {result.stderr}"
            
            # Compare which specimens were matched in both runs
            original_specimens = set()
            rc_specimens = set()
            
            for output_dir, specimen_set in [(original_output, original_specimens), (rc_output, rc_specimens)]:
                full_dir = output_dir / "full"
                if full_dir.exists():
                    for pool_dir in full_dir.iterdir():
                        if pool_dir.is_dir():
                            for specimen_file in pool_dir.glob("*.fastq"):
                                # Extract specimen name from filename (no prefix to remove)
                                specimen_name = specimen_file.stem
                                specimen_set.add(specimen_name)
            
            # Both runs should identify the same specimens (same biological matches)
            common_specimens = original_specimens.intersection(rc_specimens)
            
            assert len(common_specimens) > 0, (
                f"No common specimens found between orientations.\n"
                f"Original: {original_specimens}\n"
                f"RC: {rc_specimens}\n"
                f"This suggests the matching logic may not be working correctly in both orientations."
            )
            
            print(f"Successfully matched {len(common_specimens)} specimens in both orientations: {sorted(common_specimens)}")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])