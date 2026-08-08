"""Property tests for the prefix inverted index prefilter.

The index must be a strict superset of what SHW (prefix) alignment can
accept: any barcode whose reverse complement aligns to a window within
max_distance MUST appear in candidates() for that window. A false
negative here silently drops reads, so this is the load-bearing property.
"""

import random

import edlib
from Bio.Seq import reverse_complement

from specimux.prefix_index import PrefixBarcodePrefilter

ALPHABET = "ACGT"


def _random_barcodes(n, length, rng):
    seen = set()
    while len(seen) < n:
        seen.add("".join(rng.choice(ALPHABET) for _ in range(length)))
    return sorted(seen)


def _mutate(s, n_edits, rng):
    s = list(s)
    for _ in range(n_edits):
        op = rng.choice("sid")
        pos = rng.randrange(len(s))
        if op == "s":
            s[pos] = rng.choice(ALPHABET)
        elif op == "i":
            s.insert(pos, rng.choice(ALPHABET))
        elif op == "d" and len(s) > 1:
            del s[pos]
    return "".join(s)


def test_no_false_negatives():
    rng = random.Random(1234)
    max_dist = 3
    barcodes = _random_barcodes(24, 13, rng)
    pairs = [(b, reverse_complement(b)) for b in barcodes]
    prefilter = PrefixBarcodePrefilter.build(pairs, max_dist)

    checked = 0
    for _ in range(3000):
        b, b_rc = rng.choice(pairs)
        window = _mutate(b_rc, rng.randint(0, max_dist), rng) + \
            "".join(rng.choice(ALPHABET) for _ in range(30))
        r = edlib.align(b_rc, window, "SHW", "distance", max_dist)
        if r["editDistance"] == -1:
            continue
        checked += 1
        cands = prefilter.candidates(window, 0)
        assert (b, b_rc) in cands, (
            f"false negative: {b_rc} SHW-matches window at distance "
            f"{r['editDistance']} but is not in candidates")
    assert checked > 1000  # the test must actually exercise matches


def test_no_false_negatives_all_barcodes_vs_window():
    """Every barcode that aligns to a shared window must be a candidate."""
    rng = random.Random(99)
    max_dist = 3
    barcodes = _random_barcodes(40, 13, rng)
    pairs = [(b, reverse_complement(b)) for b in barcodes]
    prefilter = PrefixBarcodePrefilter.build(pairs, max_dist)

    for _ in range(500):
        window = "".join(rng.choice(ALPHABET) for _ in range(40))
        cands = set(prefilter.candidates(window, 0))
        for b, b_rc in pairs:
            r = edlib.align(b_rc, window, "SHW", "distance", max_dist)
            if r["editDistance"] != -1:
                assert (b, b_rc) in cands


def test_deterministic_build():
    rng = random.Random(7)
    barcodes = _random_barcodes(12, 13, rng)
    pairs = [(b, reverse_complement(b)) for b in barcodes]
    a = PrefixBarcodePrefilter.build(pairs, 3)
    b = PrefixBarcodePrefilter.build(pairs, 3)
    assert a._index == b._index


def test_mixed_length_barcodes_degrade_gracefully():
    pairs = [("ACGTACGTACGTA", reverse_complement("ACGTACGTACGTA")),
             ("ACGTACGT", reverse_complement("ACGTACGT"))]
    assert PrefixBarcodePrefilter.try_load_or_build(pairs, 3) is None
