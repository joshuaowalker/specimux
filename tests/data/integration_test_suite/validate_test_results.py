#!/usr/bin/env python3
"""
Validation script for specimux integration test results.
Compares actual output against expected output.

This validation is order-independent and suitable for multiprocessing environments.
It validates directory structure, file counts, and sequence assignments without 
depending on file ordering within directories.
"""

import sys
import os
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

def count_files_in_directory(directory):
    """Count files recursively in directory structure."""
    if not os.path.exists(directory):
        return {}
    
    counts = defaultdict(int)
    for root, dirs, files in os.walk(directory):
        # Get relative path for comparison
        rel_path = os.path.relpath(root, directory)
        if rel_path == '.':
            rel_path = 'root'
        
        # Count .fastq files only
        fastq_files = [f for f in files if f.endswith('.fastq')]
        if fastq_files:
            counts[rel_path] = len(fastq_files)
    
    return dict(counts)

def get_sequence_assignments(directory):
    """Get mapping of sequence IDs to their file paths with counts.
    
    Returns a dictionary where each sequence ID maps to a dictionary of
    {file_path: count} showing how many times the sequence appears in each file.
    """
    if not os.path.exists(directory):
        return {}
    
    assignments = defaultdict(lambda: defaultdict(int))
    for root, dirs, files in os.walk(directory):
        # Skip trace directory
        if 'trace' in root:
            continue
            
        for file in files:
            if file.endswith('.fastq'):
                file_path = os.path.join(root, file)
                rel_path = os.path.relpath(file_path, directory)
                
                # Parse FASTQ and collect sequence IDs
                try:
                    for record in SeqIO.parse(file_path, "fastq"):
                        assignments[record.id][rel_path] += 1
                except Exception as e:
                    print(f"Warning: Could not parse {file_path}: {e}")
    
    return dict(assignments)

def compare_output_structure(actual_dir, expected_dir):
    """Compare the structure and file counts between actual and expected output."""
    
    print(f"Comparing output structure:")
    print(f"  Actual: {actual_dir}")
    print(f"  Expected: {expected_dir}")
    print()
    
    actual_counts = count_files_in_directory(actual_dir)
    expected_counts = count_files_in_directory(expected_dir)
    
    # Check for matching structure
    all_paths = set(actual_counts.keys()) | set(expected_counts.keys())
    
    success = True
    for path in sorted(all_paths):
        actual_count = actual_counts.get(path, 0)
        expected_count = expected_counts.get(path, 0)
        
        status = "✓" if actual_count == expected_count else "✗"
        if actual_count != expected_count:
            success = False
            
        print(f"{status} {path:40} | Actual: {actual_count:3d} | Expected: {expected_count:3d}")
    
    return success

def compare_sequence_assignments(actual_dir, expected_dir):
    """Compare sequence assignments between actual and expected output.
    
    Handles sequences that appear in multiple locations by comparing the
    complete set of locations and counts for each sequence.
    """
    print(f"\nComparing sequence assignments:")
    
    actual_assignments = get_sequence_assignments(actual_dir)
    expected_assignments = get_sequence_assignments(expected_dir)
    
    # Get all sequence IDs from both sets
    all_sequences = set(actual_assignments.keys()) | set(expected_assignments.keys())
    
    correct_assignments = 0
    total_assignments = len(all_sequences)
    
    mismatches = []
    
    for seq_id in sorted(all_sequences):
        actual_locations = actual_assignments.get(seq_id, {})
        expected_locations = expected_assignments.get(seq_id, {})
        
        # Compare the complete set of locations and counts
        if actual_locations == expected_locations:
            correct_assignments += 1
        else:
            mismatches.append((seq_id, actual_locations, expected_locations))
    
    print(f"  Sequence assignments: {correct_assignments}/{total_assignments} correct")
    
    if mismatches:
        print(f"  Mismatched assignments ({len(mismatches)}):")
        for seq_id, actual_locs, expected_locs in mismatches[:10]:  # Show first 10
            print(f"    {seq_id}")
            
            # Format actual locations
            if not actual_locs:
                print(f"      Actual:   MISSING")
            elif len(actual_locs) == 1:
                file_path, count = list(actual_locs.items())[0]
                print(f"      Actual:   {file_path}" + (f" (x{count})" if count > 1 else ""))
            else:
                print(f"      Actual:   {len(actual_locs)} locations:")
                for file_path, count in sorted(actual_locs.items()):
                    print(f"                {file_path}" + (f" (x{count})" if count > 1 else ""))
            
            # Format expected locations
            if not expected_locs:
                print(f"      Expected: MISSING")
            elif len(expected_locs) == 1:
                file_path, count = list(expected_locs.items())[0]
                print(f"      Expected: {file_path}" + (f" (x{count})" if count > 1 else ""))
            else:
                print(f"      Expected: {len(expected_locs)} locations:")
                for file_path, count in sorted(expected_locs.items()):
                    print(f"                {file_path}" + (f" (x{count})" if count > 1 else ""))
                    
        if len(mismatches) > 10:
            print(f"    ... and {len(mismatches) - 10} more")
    
    # Report statistics on multiple assignments
    actual_multiples = sum(1 for locs in actual_assignments.values() if len(locs) > 1)
    expected_multiples = sum(1 for locs in expected_assignments.values() if len(locs) > 1)
    
    if actual_multiples > 0 or expected_multiples > 0:
        print(f"\n  Sequences with multiple locations:")
        print(f"    Actual:   {actual_multiples} sequences")
        print(f"    Expected: {expected_multiples} sequences")
    
    return correct_assignments == total_assignments

def check_trace_files(actual_dir, expected_dir):
    """Check that trace files were generated."""
    actual_trace = Path(actual_dir) / "trace"
    expected_trace = Path(expected_dir) / "trace"
    
    actual_trace_files = list(actual_trace.glob("*.tsv")) if actual_trace.exists() else []
    expected_trace_files = list(expected_trace.glob("*.tsv")) if expected_trace.exists() else []
    
    print(f"\nTrace files:")
    print(f"  Actual trace files: {len(actual_trace_files)}")
    print(f"  Expected trace files: {len(expected_trace_files)}")
    
    return len(actual_trace_files) > 0 and len(expected_trace_files) > 0

def main():
    if len(sys.argv) != 3:
        print("Usage: python validate_test_results.py <actual_output_dir> <expected_output_dir>")
        sys.exit(1)
    
    actual_dir = sys.argv[1]
    expected_dir = sys.argv[2]
    
    if not os.path.exists(actual_dir):
        print(f"Error: Actual output directory does not exist: {actual_dir}")
        sys.exit(1)
        
    if not os.path.exists(expected_dir):
        print(f"Error: Expected output directory does not exist: {expected_dir}")
        sys.exit(1)
    
    print("Specimux Integration Test Validation")
    print("=" * 50)
    
    # Compare output structure
    structure_ok = compare_output_structure(actual_dir, expected_dir)
    
    # Compare sequence assignments
    assignments_ok = compare_sequence_assignments(actual_dir, expected_dir)
    
    # Check trace files
    trace_ok = check_trace_files(actual_dir, expected_dir)
    
    print("\n" + "=" * 50)
    if structure_ok and assignments_ok and trace_ok:
        print("✓ All tests PASSED")
        print("  - Output structure matches expected")
        print("  - Sequence assignments are correct")
        print("  - Trace files generated successfully")
        sys.exit(0)
    else:
        print("✗ Some tests FAILED")
        if not structure_ok:
            print("  - Output structure does not match expected")
        if not assignments_ok:
            print("  - Sequence assignments are incorrect")
        if not trace_ok:
            print("  - Trace files missing or not generated")
        sys.exit(1)

if __name__ == "__main__":
    main()