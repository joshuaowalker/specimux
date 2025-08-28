#!/usr/bin/env python3

"""
Specimen File Converter

Converts old specimen file format (with primer sequences) to new format with:
1. A primers.fasta file containing unique primer sequences
2. A new specimen file with PrimerPool column and primer names instead of sequences

Author: Josh Walker
Date: March 8, 2025
"""

import argparse
import csv
import os
from collections import defaultdict


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Convert old specimen file format to new format with primers.fasta")

    parser.add_argument("input_file", help="Input specimen file (old format)")
    parser.add_argument("--output-specimen", help="Output specimen file (new format)", default="specimen_new.txt")
    parser.add_argument("--output-primers", help="Output primers fasta file", default="primers.fasta")
    parser.add_argument("--pool-name", help="Name for the primer pool", default="pool1")

    return parser.parse_args()


def extract_primers_from_specimen_file(input_file):
    """
    Extract unique forward and reverse primers from the specimen file.

    Returns:
        tuple: (forward_primers, reverse_primers) where each is a set of unique sequences
    """
    forward_primers = set()
    reverse_primers = set()

    with open(input_file, 'r', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Skip wildcards like '*' or '-'
            if row['FwPrimer'] not in ['*', '-']:
                forward_primers.add(row['FwPrimer'])
            if row['RvPrimer'] not in ['*', '-']:
                reverse_primers.add(row['RvPrimer'])

    return forward_primers, reverse_primers


def create_primers_fasta(forward_primers, reverse_primers, output_file, pool_name):
    """
    Create a primers.fasta file with generic names for unique primer sequences.

    Returns:
        tuple: (fw_primer_map, rev_primer_map) dictionaries mapping sequences to names
    """
    fw_primer_map = {}  # Maps sequence to name
    rev_primer_map = {}  # Maps sequence to name

    fw_idx = 1
    rev_idx = 1

    with open(output_file, 'w') as f:
        # Write forward primers
        for primer_seq in sorted(forward_primers):
            primer_name = f"primer_F{fw_idx}"
            fw_primer_map[primer_seq] = primer_name
            f.write(f">{primer_name} position=forward pool={pool_name}\n")
            f.write(f"{primer_seq}\n")
            fw_idx += 1

        # Write reverse primers
        for primer_seq in sorted(reverse_primers):
            primer_name = f"primer_R{rev_idx}"
            rev_primer_map[primer_seq] = primer_name
            f.write(f">{primer_name} position=reverse pool={pool_name}\n")
            f.write(f"{primer_seq}\n")
            rev_idx += 1

    print(f"Created primers.fasta with {fw_idx - 1} forward and {rev_idx - 1} reverse primers")
    return fw_primer_map, rev_primer_map


def create_new_specimen_file(input_file, output_file, fw_primer_map, rev_primer_map, pool_name):
    """Create a new specimen file with PrimerPool column and primer names."""
    with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')

        # Create new fieldnames with PrimerPool added after SampleID
        old_fieldnames = reader.fieldnames
        new_fieldnames = ['SampleID', 'PrimerPool'] + old_fieldnames[1:]

        writer = csv.DictWriter(outfile, fieldnames=new_fieldnames, delimiter='\t')
        writer.writeheader()

        for row in reader:
            new_row = {'SampleID': row['SampleID'], 'PrimerPool': pool_name}

            # Copy existing fields except primers
            for field in old_fieldnames:
                if field != 'SampleID':
                    new_row[field] = row[field]

            # Replace primer sequences with names, preserving wildcards
            if row['FwPrimer'] in ['*', '-']:
                new_row['FwPrimer'] = row['FwPrimer']
            else:
                new_row['FwPrimer'] = fw_primer_map.get(row['FwPrimer'], row['FwPrimer'])

            if row['RvPrimer'] in ['*', '-']:
                new_row['RvPrimer'] = row['RvPrimer']
            else:
                new_row['RvPrimer'] = rev_primer_map.get(row['RvPrimer'], row['RvPrimer'])

            writer.writerow(new_row)

    print(f"Created new specimen file: {output_file}")


def main():
    args = parse_arguments()

    # Extract unique primers
    forward_primers, reverse_primers = extract_primers_from_specimen_file(args.input_file)

    # Create primers.fasta and get primer name mappings
    fw_primer_map, rev_primer_map = create_primers_fasta(
        forward_primers, reverse_primers, args.output_primers, args.pool_name)

    # Create new specimen file
    create_new_specimen_file(
        args.input_file, args.output_specimen, fw_primer_map, rev_primer_map, args.pool_name)


if __name__ == "__main__":
    main()