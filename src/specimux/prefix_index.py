#!/usr/bin/env python3

"""Prefix inverted index for barcode prefiltering.

Replaces the bloom-filter prefilter with an exact inverted index:
for each barcode we enumerate every (L-d)-mer P that could begin a
window which SHW-aligns to the barcode's reverse complement within
edit distance d, i.e. min_j edit(P, b_rc[:j]) <= d. At query time a
single dict lookup on the window prefix returns the only barcodes
worth aligning.

The index is a strict superset of what prefix alignment can accept
(restriction of an optimal alignment path to the first L-d window
characters costs at most the full-alignment distance), so it can
produce no false negatives.
"""

import hashlib
import logging
import os
import pickle
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_ALPHABET = "ACGT"
_CACHE_FORMAT_VERSION = 1


def barcode_pairs_for_prefilter(specimens) -> List[Tuple[str, str]]:
    """Collect all distinct (barcode, barcode_rc) pairs in canonical order."""
    from .constants import Primer

    seen = set()
    pairs = []
    for direction in (Primer.FWD, Primer.REV):
        for primer in specimens.get_primers(direction):
            for b, b_rc in primer.barcode_pairs:
                if b not in seen:
                    seen.add(b)
                    pairs.append((b, b_rc))
    pairs.sort(key=lambda p: p[0])
    return pairs


def _enumerate_prefixes(target: str, prefix_len: int, max_dist: int) -> List[str]:
    """All prefix_len-mers P with min_j edit(P, target[:j]) <= max_dist.

    DFS over the 4-ary trie of prefixes carrying the Levenshtein DP row
    against target; a branch is pruned as soon as min(row) exceeds
    max_dist, which is sound because extending P can never decrease the
    row minimum.
    """
    n = len(target)
    results = []
    stack = [("", list(range(n + 1)))]
    while stack:
        prefix, row = stack.pop()
        if len(prefix) == prefix_len:
            results.append(prefix)
            continue
        for c in _ALPHABET:
            prev = row[0]
            new_row = [prev + 1]
            append = new_row.append
            best = prev + 1
            for j in range(1, n + 1):
                cost = prev if target[j - 1] == c else prev + 1
                prev = row[j]
                cell = min(new_row[j - 1] + 1, prev + 1, cost)
                append(cell)
                if cell < best:
                    best = cell
            if best <= max_dist:
                stack.append((prefix + c, new_row))
    return results


class PrefixBarcodePrefilter:
    """Inverted index from window prefixes to candidate barcodes."""

    def __init__(self, barcode_pairs: List[Tuple[str, str]], max_distance: int,
                 index: Dict[str, tuple]):
        self.barcode_pairs = barcode_pairs
        self.max_distance = max_distance
        self.prefix_len = len(barcode_pairs[0][0]) - max_distance
        self._index = index

    def candidates(self, sequence: str, start: int) -> tuple:
        """Barcodes worth aligning at a window beginning at start."""
        return self._index.get(sequence[start:start + self.prefix_len], ())

    @staticmethod
    def build(barcode_pairs: List[Tuple[str, str]], max_distance: int) -> 'PrefixBarcodePrefilter':
        lengths = {len(b) for b, _ in barcode_pairs}
        if len(lengths) != 1:
            raise ValueError(f"Prefix index requires uniform barcode length, got {sorted(lengths)}")
        prefix_len = lengths.pop() - max_distance
        if prefix_len < 1:
            raise ValueError("max_distance too large for barcode length")

        t0 = time.time()
        index: Dict[str, list] = {}
        for b, b_rc in barcode_pairs:
            for p in _enumerate_prefixes(b_rc, prefix_len, max_distance):
                index.setdefault(p, []).append((b, b_rc))
        frozen = {p: tuple(v) for p, v in index.items()}
        logging.info(f"Built prefix index: {len(barcode_pairs)} barcodes, "
                     f"{len(frozen)} keys in {time.time() - t0:.1f}s")
        return PrefixBarcodePrefilter(barcode_pairs, max_distance, frozen)

    @staticmethod
    def get_cache_path(barcode_pairs: List[Tuple[str, str]], max_distance: int) -> Path:
        digest = hashlib.sha256()
        digest.update(f"v{_CACHE_FORMAT_VERSION}:d{max_distance}:".encode())
        for b, _ in barcode_pairs:
            digest.update(b.encode())
            digest.update(b",")
        cache_dir = Path.home() / ".specimux" / "cache"
        return cache_dir / f"prefix_index_{digest.hexdigest()[:24]}.pkl"

    @classmethod
    def try_load_or_build(cls, barcode_pairs: List[Tuple[str, str]],
                          max_distance: int) -> 'Optional[PrefixBarcodePrefilter]':
        """load_or_build, degrading to None (no prefilter) on ineligible input."""
        try:
            return cls.load_or_build(barcode_pairs, max_distance)
        except ValueError as e:
            logging.warning(f"Prefix index prefilter unavailable ({e}); "
                            f"barcode matching will run without prefiltering")
            return None

    @classmethod
    def load_or_build(cls, barcode_pairs: List[Tuple[str, str]],
                      max_distance: int) -> 'PrefixBarcodePrefilter':
        """Load the index from cache, building and caching it on a miss."""
        cache_path = cls.get_cache_path(barcode_pairs, max_distance)
        if cache_path.exists():
            try:
                with open(cache_path, "rb") as fh:
                    index = pickle.load(fh)
                logging.info(f"Loaded prefix index cache: {cache_path}")
                return cls(barcode_pairs, max_distance, index)
            except (OSError, pickle.UnpicklingError) as e:
                logging.warning(f"Could not load prefix index cache ({e}); rebuilding")

        prefilter = cls.build(barcode_pairs, max_distance)
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            fd, tmp = tempfile.mkstemp(dir=cache_path.parent, suffix=".tmp")
            with os.fdopen(fd, "wb") as fh:
                pickle.dump(prefilter._index, fh, protocol=pickle.HIGHEST_PROTOCOL)
            os.replace(tmp, cache_path)
            logging.info(f"Cached prefix index: {cache_path}")
        except OSError as e:
            logging.warning(f"Could not cache prefix index: {e}")
        return prefilter
