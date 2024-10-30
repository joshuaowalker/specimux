#!/usr/bin/env python3.10

"""
specimine.py: A tool for mining additional matches from raw reads using reference sequences.

This tool searches for specific sequences in FASTQ/FASTA files, excluding certain sequences,
and outputs matched sequences to separate files. It can operate in manual or auto mode.

Usage (manual mode):
    python specimine.py input.fastq search.fasta -e exclude.fasta -i 0.9 -o output_dir

Usage (auto mode):
    python specimine.py --auto --index-file Index.txt --summary-file summary.fasta [--forward] [--reverse] [--max-ric 100]

Dependencies:
    - Biopython
    - edlib

To install dependencies:
    pip install biopython edlib
"""

import argparse
import logging
import os
import re
import sys
from re import search
from typing import List, Dict, Tuple

import edlib
from Bio import SeqIO
from Bio.Seq import reverse_complement

def parse_arguments():
    parser = argparse.ArgumentParser(description="Mine additional specimen matches from raw reads.")
    parser.add_argument("--auto", action="store_true", help="Enable auto mode")
    parser.add_argument("--index-file", help="Path to the Index.txt file (required for auto mode)")
    parser.add_argument("--summary-file", help="Path to the summary.fasta file (required for auto mode)")
    parser.add_argument("--forward", action="store_true", help="Enable forward barcodes in auto mode")
    parser.add_argument("--reverse", action="store_true", default=True, help="Enable reverse barcodes in auto mode")
    parser.add_argument("--max-ric", type=int, default=100, help="Maximum RIC value for search sequences in auto mode")
    parser.add_argument("--write-augmented", action="store_true",
                        help="Write augmented FASTQ files containing original and mined sequences")
    parser.add_argument("input_file", nargs="?", help="Input FASTQ/FASTA file to be mined (for manual mode)")
    parser.add_argument("search_file", nargs="?", help="FASTA file of sequences to search for (for manual mode)")
    parser.add_argument("-e", "--exclude_file", help="Optional FASTA file of sequences to exclude (for manual mode)")
    parser.add_argument("-i", "--min_identity", type=float, default=0.9,
                        help="Minimum alignment identity for a match (default: 0.9)")
    parser.add_argument("-o", "--output_dir", default=".", help="Output directory for mined sequences")

    parser.add_argument("--debug", action="store_true", help="Enable debug mode")
    return parser.parse_args()


def read_fasta(filename: str) -> Dict[str, str]:
    return {record.id: str(record.seq) for record in SeqIO.parse(filename, "fasta")}


def calculate_identity(alignment_result: Dict, query_length: int) -> float:
    d = alignment_result["editDistance"]
    logging.debug(f"aligned with {d} edits")
    if d == -1:
        return 0
    return 1 - (d / query_length)


def check_sequence_similarities(sequences: Dict[str, str], min_identity: float) -> List[Tuple[str, str, float]]:
    similarities = []
    for id1, seq1 in sequences.items():
        for id2, seq2 in sequences.items():
            if id1 < id2:
                max_distance = int(len(seq1) * (1 - min_identity))
                alignment = edlib.align(seq1, seq2, mode="HW", task="path", k=max_distance)
                identity = calculate_identity(alignment, len(seq1))
                if identity >= min_identity:
                    similarities.append((id1, id2, identity))
    return similarities


def preprocess_sequences(search_seqs: Dict[str, str], exclude_seqs: Dict[str, str], min_identity: float) -> Tuple[
    Dict[str, str], Dict[str, str]]:
    # Add primers to search and exclude sequences
    search_seqs = {id: FWD_PRIMER + seq + REV_PRIMER for id, seq in search_seqs.items()}
    exclude_seqs = {id: FWD_PRIMER + seq + REV_PRIMER for id, seq in exclude_seqs.items()}

    # Check search sequences against each other
    search_similarities = check_sequence_similarities(search_seqs, min_identity)
    for id1, id2, identity in search_similarities:
        logging.warning(f"Search sequences {id1} and {id2} align with identity {identity:.2f}")

    # Check search sequences against exclusion sequences
    if exclude_seqs:
        unexclude = set()
        for search_id, search_seq in search_seqs.items():
            for exclude_id, exclude_seq in exclude_seqs.items():
                max_distance = int(len(search_seq) * (1 - min_identity))
                alignment = edlib.align(search_seq, exclude_seq, mode="HW", task="path", k=max_distance)
                identity = calculate_identity(alignment, len(search_seq))
                if identity >= min_identity:
                    if identity == 1.0:
                        logging.info(f"Removing 100% matching exclusion sequence {exclude_id}")
                        unexclude.add(exclude_id)
                    else:
                        logging.warning(
                            f"Search sequence {search_id} aligns with exclusion sequence {exclude_id} (identity: {identity:.2f})")

        for i in unexclude:
            del exclude_seqs[i]
    return search_seqs, exclude_seqs


def process_sequences(input_file: str, search_seqs: Dict[str, str], exclude_seqs: Dict[str, str], min_identity: float,
                      output_dir: str, write_augmented: bool = False):
    os.makedirs(output_dir, exist_ok=True)
    input_format = "fastq" if input_file.endswith((".fastq", ".fq")) else "fasta"
    output_files = {id: open(os.path.join(output_dir, f"sample_{id}_mined.{input_format}"), "w") for id in search_seqs}

    if write_augmented:
        # Group output files by specimen ID (everything before first "-")
        specimen_files = {}
        for search_id in search_seqs:
            specimen_id = re.match(r'(.*)-\d+', search_id).group(1)
            if specimen_id not in specimen_files:
                specimen_files[specimen_id] = set()
            specimen_files[specimen_id].add(search_id)

        # Create augmented files
        augmented_files = {}
        for specimen_id in specimen_files:
            augmented_filename = os.path.join(output_dir, f"sample_{specimen_id}_augmented.{input_format}")
            augmented_files[specimen_id] = open(augmented_filename, "w")

            # Copy original file contents to augmented file
            original_filename = f"sample_{specimen_id}.{input_format}"
            if os.path.exists(original_filename):
                with open(original_filename) as f:
                    augmented_files[specimen_id].write(f.read())
            else:
                logging.error(f"Could not find {original_filename} for --write-augmented")

    total_reads = 0
    matched_reads = 0
    excluded_reads = 0

    for record in SeqIO.parse(input_file, input_format):
        total_reads += 1
        seq = str(record.seq)

        # Check against exclusion sequences
        exclude_match = False
        for exclude_id, exclude_seq in exclude_seqs.items():
            max_distance = int(len(exclude_seq) * (1 - min_identity))
            alignment = edlib.align(exclude_seq, seq, mode="HW", task="path", k=max_distance)
            if alignment["editDistance"] != -1 and calculate_identity(alignment, len(exclude_seq)) >= min_identity:
                exclude_match = True
                excluded_reads += 1
                break

        if exclude_match:
            continue

        # Check against search sequences
        best_match = None
        best_identity = 0
        best_alignment = None
        for search_id, search_seq in search_seqs.items():
            max_distance = int(len(search_seq) * (1 - min_identity))
            alignment = edlib.align(search_seq, seq, mode="HW", task="path", k=max_distance)
            if alignment["editDistance"] != -1:
                identity = calculate_identity(alignment, len(search_seq))
                if identity >= min_identity and identity > best_identity:
                    best_match = search_id
                    best_identity = identity
                    best_alignment = alignment

        if best_match:
            # Trim the sequence to the matched region
            start, end = best_alignment["locations"][0]
            trimmed_seq = seq[start:end + 1]
            trimmed_record = record[start:end + 1]

            # Write the trimmed sequence to mined file
            SeqIO.write(trimmed_record, output_files[best_match], input_format)

            # Also write to augmented file if enabled
            if write_augmented:
                specimen_id = re.match(r'(.*)-\d+', search_id).group(1)
                SeqIO.write(trimmed_record, augmented_files[specimen_id], input_format)

            matched_reads += 1

        # Update progress
        sys.stdout.write(f"\rProcessed: {total_reads}, Matched: {matched_reads}, Excluded: {excluded_reads}")
        sys.stdout.flush()

    sys.stdout.write("\n")  # Move to the next line after progress indicator

    # Close all files
    for file in output_files.values():
        file.close()
    if write_augmented:
        for file in augmented_files.values():
            file.close()


FWD_PRIMER = "CTTGGTCATTTAGAGGAAGTAA"
REV_PRIMER = reverse_complement("TCCTCCGCTTATTGATATGC")  # GCATATCAATAAGCGGAGGA

def parse_index_file(index_file: str, forward, reverse) -> Dict[str, Dict[str, str]]:
    samples = {}
    with open(index_file, 'r') as f:
        next(f)  # Skip header
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                sample_id, forward_barcode, _, reverse_barcode, _ = fields[:5]
                if forward:
                    key = ("fwd", forward_barcode)
                    if key not in samples:
                        samples[key] = []
                    samples[key].append(sample_id)
                if reverse:
                    key = ("rev", reverse_barcode)
                    if key not in samples:
                        samples[key] = []
                    samples[key].append(sample_id)
    return samples

def parse_summary_fasta(summary_file: str) -> Dict[str, Tuple[str, int]]:
    sequences = {}
    for record in SeqIO.parse(summary_file, "fasta"):
        ric_match = re.search(r'ric=(\d+)', record.description)
        if ric_match:
            ric = int(ric_match.group(1))
            sequences[record.id] = (str(record.seq), ric)
    return sequences

def auto_mode(args):
    samples = parse_index_file(args.index_file, args.forward, args.reverse)
    all_sequences = parse_summary_fasta(args.summary_file)

    for key, samples in samples.items():
        (which, barcode) = key
        if which == "fwd" and args.forward:
            process_barcode_file(f"sample_barcode_fwd_{barcode}.fastq", samples, all_sequences, args)
        elif which == "rev" and args.reverse:
            process_barcode_file(f"sample_barcode_rev_{barcode}.fastq", samples, all_sequences, args)


def process_barcode_file(input_file: str, samples: List[str], all_sequences: Dict[str, Tuple[str, int]], args):
    search_seqs = {}
    exclude_seqs = {}

    for seq_id, (sequence, ric) in all_sequences.items():
        for sample_id in samples:
            if seq_id.startswith(sample_id):
                if ric < args.max_ric:
                    search_seqs[seq_id] = sequence
                else:
                    exclude_seqs[seq_id] = sequence

    print(
        f"Processing {input_file} for {len(samples)} samples: {len(search_seqs)} search sequences, {len(exclude_seqs)} exclude sequences")
    if len(search_seqs) > 0:
        search_seqs, exclude_seqs = preprocess_sequences(search_seqs, exclude_seqs, args.min_identity)
        process_sequences(input_file, search_seqs, exclude_seqs, args.min_identity, args.output_dir,
                          args.write_augmented)

def main():
    args = parse_arguments()
    l = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=l, format='%(asctime)s - %(levelname)s - %(message)s')

    if args.auto:
        if not (args.index_file and args.summary_file):
            logging.error("In auto mode, both --index-file and --summary-file are required.")
            sys.exit(1)
        auto_mode(args)
    else:
        if not (args.input_file and args.search_file):
            logging.error("In manual mode, input_file and search_file are required.")
            sys.exit(1)

        search_seqs = read_fasta(args.search_file)
        exclude_seqs = read_fasta(args.exclude_file) if args.exclude_file else {}

        search_seqs, exclude_seqs = preprocess_sequences(search_seqs, exclude_seqs, args.min_identity)
        process_sequences(args.input_file, search_seqs, exclude_seqs, args.min_identity, args.output_dir)


if __name__ == "__main__":
    main()


