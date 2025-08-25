#!/usr/bin/env python3
"""
Integration test runner for specimux.
Runs the integration test suite and validates results.
"""

import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

def run_integration_test():
    """Run specimux integration test and validate results."""
    
    # Test data paths
    test_dir = Path(__file__).parent / "data" / "integration_test_suite"
    primers_file = test_dir / "primers.fasta"
    specimens_file = test_dir / "specimens.txt"
    sequences_file = test_dir / "sequences.fastq"
    expected_output = test_dir / "expected_output"
    validation_script = test_dir / "validate_test_results.py"
    
    # Check test files exist
    required_files = [primers_file, specimens_file, sequences_file, expected_output, validation_script]
    for file_path in required_files:
        if not file_path.exists():
            print(f"Error: Required test file missing: {file_path}")
            return False
    
    # Create temporary output directory
    with tempfile.TemporaryDirectory(prefix="specimux_test_") as temp_dir:
        output_dir = Path(temp_dir) / "output"
        
        print("Running specimux integration test...")
        print(f"Input: {sequences_file} ({40} sequences)")
        print(f"Output: {output_dir}")
        print()
        
        # Run specimux
        cmd = [
            sys.executable, "-m", "specimux.cli",
            str(primers_file),
            str(specimens_file), 
            str(sequences_file),
            "-F", "-O", str(output_dir), "-d"
        ]
        
        try:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print("✓ Specimux completed successfully")
        except subprocess.CalledProcessError as e:
            print("✗ Specimux failed:")
            print(f"  Command: {' '.join(cmd)}")
            print(f"  Return code: {e.returncode}")
            print(f"  Error output: {e.stderr}")
            return False
        
        # Run validation
        print("\nValidating results...")
        validation_cmd = [
            sys.executable, str(validation_script),
            str(output_dir),
            str(expected_output)
        ]
        
        try:
            validation_result = subprocess.run(validation_cmd, check=True, capture_output=True, text=True)
            print(validation_result.stdout)
            return True
        except subprocess.CalledProcessError as e:
            print("✗ Validation failed:")
            print(e.stdout)
            return False

def main():
    """Main entry point."""
    print("Specimux Integration Test Runner")
    print("=" * 40)
    print()
    
    
    success = run_integration_test()
    
    print("\n" + "=" * 40)
    if success:
        print("🎉 Integration tests PASSED")
        sys.exit(0)
    else:
        print("💥 Integration tests FAILED")
        sys.exit(1)

if __name__ == "__main__":
    main()