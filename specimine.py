#!/usr/bin/env python3.10

"""
specimine.py: A tool for mining additional matches from raw reads using reference sequences.

This tool searches for specific sequences in FASTQ/FASTA files, excluding certain sequences,
and outputs matched sequences to separate files.

Usage:
    python specimine.py input.fastq search.fasta -e exclude.fasta -i 0.9 -o output_dir

Dependencies:
    - Biopython
    - edlib

To install dependencies:
    pip install biopython edlib
"""

import argparse
import os
import sys
import logging
from typing import List, Dict, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
import edlib
from Bio.Seq import reverse_complement


def parse_arguments():
    parser = argparse.ArgumentParser(description="Mine fungal ITS amplicons from raw reads.")
    parser.add_argument("input_file", help="Input FASTQ/FASTA file to be mined")
    parser.add_argument("search_file", help="FASTA file of sequences to search for")
    parser.add_argument("-e", "--exclude_file", help="Optional FASTA file of sequences to exclude")
    parser.add_argument("-i", "--min_identity", type=float, default=0.9,
                        help="Minimum alignment identity for a match (default: 0.9)")
    parser.add_argument("-o", "--output_dir", default=".", help="Output directory for mined sequences")
    return parser.parse_args()


def read_fasta(filename: str) -> Dict[str, str]:
    return {record.id: str(record.seq) for record in SeqIO.parse(filename, "fasta")}


def calculate_identity(alignment_result: Dict, query_length: int) -> float:
    if alignment_result["editDistance"] == -1:
        return 0
    return 1 - (alignment_result["editDistance"] / query_length)


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
        unexclude = []
        for search_id, search_seq in search_seqs.items():
            for exclude_id, exclude_seq in exclude_seqs.items():
                max_distance = int(len(search_seq) * (1 - min_identity))
                alignment = edlib.align(search_seq, exclude_seq, mode="HW", task="path", k=max_distance)
                identity = calculate_identity(alignment, len(search_seq))
                if identity >= min_identity:
                    if identity == 1.0:
                        logging.info(f"Removing 100% matching exclusion sequence {exclude_id}")
                        unexclude.append(exclude_id)
                    else:
                        logging.warning(
                            f"Search sequence {search_id} aligns with exclusion sequence {exclude_id} (identity: {identity:.2f})")

        for i in unexclude:
            del exclude_seqs[i]
    return search_seqs, exclude_seqs


def process_sequences(input_file: str, search_seqs: Dict[str, str], exclude_seqs: Dict[str, str], min_identity: float,
                      output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    input_format = "fastq" if input_file.endswith((".fastq", ".fq")) else "fasta"
    output_files = {id: open(os.path.join(output_dir, f"sample_mined_{id}.{input_format}"), "w") for id in search_seqs}

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

            # Write the trimmed sequence
            SeqIO.write(trimmed_record, output_files[best_match], input_format)
            matched_reads += 1

        # Update progress
        sys.stdout.write(f"\rProcessed: {total_reads}, Matched: {matched_reads}, Excluded: {excluded_reads}")
        sys.stdout.flush()

    sys.stdout.write("\n")  # Move to the next line after progress indicator
    for file in output_files.values():
        file.close()


FWD_PRIMER = "CTTGGTCATTTAGAGGAAGTAA"
REV_PRIMER = reverse_complement("TCCTCCGCTTATTGATATGC")  # GCATATCAATAAGCGGAGGA


def main():
    args = parse_arguments()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    search_seqs = read_fasta(args.search_file)
    exclude_seqs = read_fasta(args.exclude_file) if args.exclude_file else {}

    search_seqs, exclude_seqs = preprocess_sequences(search_seqs, exclude_seqs, args.min_identity)
    process_sequences(args.input_file, search_seqs, exclude_seqs, args.min_identity, args.output_dir)


if __name__ == "__main__":
    main()