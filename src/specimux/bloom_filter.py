#!/usr/bin/env python3

"""
Bloom filter implementation for barcode prefiltering.

This module contains bloom filter functionality
extracted from the original specimux.py monolithic file.
"""

import fcntl
import hashlib
import logging
import os
from typing import List, Optional, Set

import pybloomfilter
from Bio.Seq import reverse_complement
from tqdm import tqdm

from .constants import Primer
from .databases import Specimens


class BloomPrefilter:
    """Fast barcode pre-filter using Bloom filters with proper file handling"""

    def __init__(self, barcodes: List[str], max_distance: int, error_rate: float = 0.05,
                 filename: Optional[str] = None):
        """Initialize bloom filter for barcode matching

        Args:
            barcodes: List of barcodes to add to filter
            max_distance: Maximum edit distance for matches
            error_rate: Acceptable false positive rate
            filename: Optional file to store the filter
        """
        if not barcodes:
            raise ValueError("Must provide at least one barcode")

        self.barcodes = list(set(barcodes))  # Deduplicate
        self.barcode_length = len(self.barcodes[0])
        self.max_distance = max_distance
        self.error_rate = error_rate
        self.min_length = self.barcode_length - max_distance

        # Calculate filter size using first barcode as sample
        sample_variants = len(self._generate_variants(self.barcodes[0]))
        total_combinations = sample_variants * len(self.barcodes)
        logging.info(f"Building Bloom filter for estimated {total_combinations} variants")

        # Create and populate filter
        if filename:
            self.bloom_filter = pybloomfilter.BloomFilter(total_combinations, error_rate, filename)
        else:
            self.bloom_filter = pybloomfilter.BloomFilter(total_combinations, error_rate)

        total_variants = 0
        # Add progress bar wrapping barcodes iteration
        for barcode in tqdm(self.barcodes, desc="Building Bloom filter", unit="barcode"):
            variants = self._generate_variants(barcode)
            total_variants += len(variants)
            for variant in variants:
                # Use fixed width format with truncation
                truncated = variant[:self.min_length]
                key = barcode + truncated
                self.bloom_filter.add(key)

        logging.info(f"Added {total_variants} actual variants to Bloom filter")

    def _generate_variants(self, barcode: str) -> Set[str]:
        """Generate all possible variants of a barcode within max_distance."""
        variants = {barcode}
        bases = "ACGT"

        def add_variants_recursive(current: str, distance_left: int, variants: Set[str]):
            if distance_left == 0:
                return

            # Substitutions
            for i in range(len(current)):
                for base in bases:
                    if base != current[i]:
                        new_str = current[:i] + base + current[i + 1:]
                        variants.add(new_str)
                        add_variants_recursive(new_str, distance_left - 1, variants)

            # Deletions
            for i in range(len(current)):
                new_str = current[:i] + current[i + 1:]
                variants.add(new_str)
                add_variants_recursive(new_str, distance_left - 1, variants)

            # Insertions
            for i in range(len(current) + 1):
                for base in bases:
                    new_str = current[:i] + base + current[i:]
                    variants.add(new_str)
                    add_variants_recursive(new_str, distance_left - 1, variants)

        add_variants_recursive(barcode, self.max_distance, variants)
        return variants

    @staticmethod
    def get_cache_path(barcodes: List[str], max_distance: int, error_rate: float = 0.05) -> str:
        """Generate a unique cache filename based on inputs."""
        m = hashlib.sha256()
        for bc in sorted(barcodes):  # Sort for consistency
            m.update(bc.encode())
        m.update(str(max_distance).encode())
        m.update(str(error_rate).encode())

        hash_prefix = m.hexdigest()[:16]
        cache_dir = os.path.join(os.path.expanduser("~"), ".specimux", "cache")
        os.makedirs(cache_dir, exist_ok=True)

        return os.path.join(cache_dir,
                            f"bloom_prefilter_{hash_prefix}_k{max_distance}_e{error_rate}.bf")

    @classmethod
    def create_filter(cls, barcode_rcs: List[str], max_distance: int, error_rate: float = 0.05) -> str:
        """Initialize bloom filter for barcode matching

        Args:
            barcode_rcs: Barcodes in the direction which they will be matched
            max_distance: Maximum edit distance for matches
            error_rate: Acceptable false positive rate

        Returns:
            str: Path to the created/cached bloom filter file
        """

        # Get cache paths
        cache_path = cls.get_cache_path(barcode_rcs, max_distance, error_rate)
        lock_path = cache_path + '.lock'

        # Create lock file if needed
        if not os.path.exists(lock_path):
            open(lock_path, 'w').close()

        # Create filter with proper locking
        with open(lock_path, 'r') as lock_file:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
            try:
                if not os.path.exists(cache_path):
                    logging.info("Creating initial Bloom filter...")
                    prefilter = cls(barcode_rcs, max_distance, error_rate, filename=cache_path)
                    prefilter.save(cache_path)
                    prefilter.close()
                    logging.info(f"Saved Bloom filter to cache: {cache_path}")
            finally:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)

        return cache_path

    @classmethod
    def load_readonly(cls, filename: str, barcodes: List[str], max_distance: int) -> 'BloomPrefilter':
        """Load BloomPrefilter from file in read-only mode"""
        instance = cls.__new__(cls)
        instance.barcodes = list(set(barcodes))
        instance.barcode_length = len(instance.barcodes[0])
        instance.max_distance = max_distance
        instance.min_length = instance.barcode_length - max_distance

        try:
            # Open with read-only mmap
            instance.bloom_filter = pybloomfilter.BloomFilter.open(filename, mode='r')
        except Exception as e:
            raise ValueError(f"Failed to load read-only filter: {e}")

        return instance

    def save(self, filename: str):
        """Save filter to file with proper syncing"""
        self.bloom_filter.sync()

    def match(self, barcode: str, sequence: str) -> bool:
        """Check if barcode matches sequence within max_distance."""
        # Truncate sequence to minimum length
        truncated = sequence[:self.min_length]

        if barcode not in self.barcodes:
            return True

        # Concatenate in same order as when building filter
        key = barcode + truncated
        return key in self.bloom_filter

    def close(self):
        """Close the filter and release resources"""
        if hasattr(self, 'bloom_filter'):
            self.bloom_filter.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


def barcodes_for_bloom_prefilter(specimens: Specimens) -> List[str]:
    """Extract barcodes for bloom filter prefiltering."""
    # Collect barcodes
    all_b1s = set()
    all_b2s = set()
    for primer in specimens.get_primers(Primer.FWD):
        all_b1s.update(primer.barcodes)
    for primer in specimens.get_primers(Primer.REV):
        all_b2s.update(primer.barcodes)
    barcode_rcs = [reverse_complement(b) for b in all_b1s] + [reverse_complement(b) for b in all_b2s]
    return barcode_rcs