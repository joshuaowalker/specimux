#!/usr/bin/env python3

"""
specimine.py: Mine additional candidate sequences from partial matches for subsequent clustering.

This tool finds additional candidate sequences from partial barcode matches by comparing
them against known good sequences from a specimen's full matches.

Usage:
    python specimine.py --index INDEX.txt --fastq SPECIMEN.fastq [--partial-forward] [--partial-reverse] [--min-identity 0.85]

Arguments:
    --index: Path to specimen index file (same as used with specimux)
    --fastq: Path to full match FASTQ file for a specimen
    --partial-forward: Include forward partial matches (default: False)
    --partial-reverse: Include reverse partial matches (default: True)
    --min-identity: Minimum alignment identity for a match (default: 0.85)
    --debug: Enable debug logging

Dependencies:
    - Biopython
    - edlib
"""

import argparse
import logging
import os
import re
import sys
from typing import Dict, List, Tuple

import edlib
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


def parse_arguments():
    parser = argparse.ArgumentParser(description="Mine additional candidate sequences from partial matches.")
    parser.add_argument("--index", required=True, help="Path to specimen index file (same as used with specimux)")
    parser.add_argument("--fastq", required=True, help="Path to full match FASTQ file for a specimen")
    parser.add_argument("--partial-forward", action="store_true", default=False,
                        help="Include forward partial matches (default: False)")
    parser.add_argument("--partial-reverse", action="store_true", default=True,
                        help="Include reverse partial matches (default: True)")
    parser.add_argument("--min-identity", type=float, default=0.85,
                        help="Minimum alignment identity for a match (default: 0.85)")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")
    return parser.parse_args()


def extract_specimen_id(fastq_path: str) -> str:
    """Extract specimen ID from the FASTQ filename."""
    filename = os.path.basename(fastq_path)
    match = re.match(r"sample_(.+)\.fastq", filename)
    if match:
        return match.group(1)
    else:
        raise ValueError(f"Could not extract specimen ID from filename: {filename}")


def find_barcodes(specimen_id: str, index_file: str) -> Tuple[str, str]:
    """Find forward and reverse barcodes for a specimen from the index file."""
    with open(index_file, 'r') as f:
        header = next(f).strip().split('\t')

        # Find indices of important columns
        sample_idx = header.index("SampleID") if "SampleID" in header else 0
        fwd_idx = header.index("FwIndex") if "FwIndex" in header else 2
        rev_idx = header.index("RvIndex") if "RvIndex" in header else 4

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > max(sample_idx, fwd_idx, rev_idx):
                if fields[sample_idx] == specimen_id:
                    return fields[fwd_idx].upper(), fields[rev_idx].upper()

    logging.error(f"Could not find specimen {specimen_id} in index file")
    return None, None


def derive_partial_match_filenames(specimen_dir: str, fwd_barcode: str, rev_barcode: str) -> Dict[str, str]:
    """Derive filenames for partial matches based on barcodes."""
    partial_files = {}

    # Directory structure: path/to/pool/primer-pair/partial/sample_barcode_*.fastq
    partial_dir = os.path.join(os.path.dirname(specimen_dir), "partial")

    if fwd_barcode:
        fwd_file = os.path.join(partial_dir, f"sample_barcode_fwd_{fwd_barcode}.fastq")
        if os.path.exists(fwd_file):
            partial_files["forward"] = fwd_file
        else:
            logging.warning(f"Forward partial match file not found: {fwd_file}")

    if rev_barcode:
        rev_file = os.path.join(partial_dir, f"sample_barcode_rev_{rev_barcode}.fastq")
        if os.path.exists(rev_file):
            partial_files["reverse"] = rev_file
        else:
            logging.warning(f"Reverse partial match file not found: {rev_file}")

    return partial_files


def calculate_identity(alignment_result: Dict, query_length: int) -> float:
    """Calculate sequence identity from alignment result."""
    d = alignment_result["editDistance"]
    if d == -1:
        return 0
    return 1 - (d / query_length)


def mine_sequences(full_match_file: str, partial_match_files: Dict[str, str], min_identity: float) -> List[SeqRecord]:
    """Mine sequences from partial matches that align well with full matches."""
    # Load full match sequences
    full_matches = list(SeqIO.parse(full_match_file, "fastq"))
    if not full_matches:
        logging.error(f"No sequences found in full match file: {full_match_file}")
        return []

    logging.info(f"Loaded {len(full_matches)} sequences from full match file")

    # Initialize list to store mined sequences
    mined_records = []

    # Process each partial match file
    for partial_type, partial_file in partial_match_files.items():
        logging.info(f"Processing {partial_type} partial matches from {partial_file}")

        # Count total sequences in partial file for tqdm
        partial_sequences = list(SeqIO.parse(partial_file, "fastq"))
        total_partials = len(partial_sequences)
        logging.info(f"Found {total_partials} sequences in partial match file")

        match_count = 0

        # Use tqdm for progress display
        for partial_record in tqdm(partial_sequences, desc=f"Mining {partial_type} matches", unit="seq"):
            partial_seq = str(partial_record.seq)

            best_match = None
            best_identity = 0

            # Compare with each full match sequence
            for full_record in full_matches:
                full_seq = str(full_record.seq)

                # Calculate maximum edit distance based on min_identity
                max_distance = int(len(full_seq) * (1 - min_identity))

                # Align sequences
                alignment = edlib.align(full_seq, partial_seq, mode="HW", task="path", k=max_distance)

                if alignment["editDistance"] != -1:
                    identity = calculate_identity(alignment, len(full_seq))

                    if identity >= min_identity and identity > best_identity:
                        best_match = full_record
                        best_identity = identity

            if best_match:
                match_count += 1
                # Modify the ID to indicate it's a mined sequence
                partial_record.id = f"{partial_record.id}_mined_{partial_type}_{best_identity:.2f}"
                partial_record.description = f"{partial_record.description} mined_{partial_type} identity={best_identity:.2f}"
                mined_records.append(partial_record)

        logging.info(f"Matched {match_count}/{total_partials} sequences from {partial_type} partial matches")

    return mined_records


def main():
    args = parse_arguments()

    # Configure logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')

    # Extract specimen ID from FASTQ filename
    specimen_id = extract_specimen_id(args.fastq)
    logging.info(f"Processing specimen: {specimen_id}")

    # Find forward and reverse barcodes from index file
    fwd_barcode, rev_barcode = find_barcodes(specimen_id, args.index)
    if not (fwd_barcode and rev_barcode):
        sys.exit(1)

    logging.info(f"Found barcodes - Forward: {fwd_barcode}, Reverse: {rev_barcode}")

    # Determine partial match filenames
    partial_match_files = derive_partial_match_filenames(os.path.dirname(args.fastq), fwd_barcode, rev_barcode)

    # Filter by user preferences
    if not args.partial_forward and "forward" in partial_match_files:
        del partial_match_files["forward"]
    if not args.partial_reverse and "reverse" in partial_match_files:
        del partial_match_files["reverse"]

    if not partial_match_files:
        logging.error("No partial match files found or selected")
        sys.exit(1)

    # Mine sequences
    mined_records = mine_sequences(args.fastq, partial_match_files, args.min_identity)
    logging.info(f"Found {len(mined_records)} mined sequences")

    # Write output file
    output_file = f"{args.fastq}.mined"
    SeqIO.write(mined_records, output_file, "fastq")
    logging.info(f"Wrote {len(mined_records)} sequences to {output_file}")


if __name__ == "__main__":
    main()