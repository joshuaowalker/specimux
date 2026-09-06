"""create_write_operation must not mutate the CandidateMatch it reads.

Dereplication can hand the same CandidateMatch to several outputs (one per
specimen when barcodes tie). Trimming used to shift the match's coordinates
in place, so every output after the first was cut from the wrong offsets.
"""

import argparse

import pytest

from specimux.constants import Barcode, Primer, ResolutionType, TrimMode
from specimux.demultiplex import create_write_operation
from specimux.models import AlignmentResult, CandidateMatch, PrimerInfo, Read

B1, P1, INSERT, P2_RC, B2_RC = "AAAA", "GGGG", "CCCCGGGG", "TTTT", "CCCC"
SEQ = B1 + P1 + INSERT + P2_RC + B2_RC


def _hit(start, end, distance=0):
    return AlignmentResult({"editDistance": distance, "locations": [(start, end)]})


def _full_match():
    read = Read("r1", "r1", SEQ, "I" * len(SEQ))
    p1 = PrimerInfo("p1", P1, Primer.FWD, ["pool"])
    p2 = PrimerInfo("p2", "AAAA", Primer.REV, ["pool"])
    m = CandidateMatch(read, 4)
    m.set_primer_match(_hit(4, 7), p1, False, Primer.FWD)
    m.set_primer_match(_hit(16, 19), p2, False, Primer.REV)
    m.add_barcode_match(_hit(0, 3), B1, False, Barcode.B1)
    m.add_barcode_match(_hit(20, 23), "CCCC", False, Barcode.B2)
    return read, m


@pytest.mark.parametrize("trim,expected", [
    (TrimMode.PRIMERS, INSERT),
    (TrimMode.BARCODES, P1 + INSERT + P2_RC),
    (TrimMode.TAILS, SEQ),
    (TrimMode.NONE, SEQ),
])
def test_repeated_outputs_are_identical(trim, expected):
    read, m = _full_match()
    args = argparse.Namespace(trim=trim)
    ops = [create_write_operation(spec, args, read, m, ResolutionType.DEREPLICATED_FULL)
           for spec in ("specA", "specB", "specC")]
    assert [op.sequence for op in ops] == [expected] * 3
    assert [op.quality_sequence for op in ops] == ["I" * len(expected)] * 3
    # Match coordinates are untouched
    assert m.get_p1_location() == (4, 7)
    assert m.get_p2_location() == (16, 19)
    assert m.get_barcode1_location() == (0, 3)
    assert m.get_barcode2_location() == (20, 23)


def test_reported_locations_are_relative_to_trimmed_sequence():
    read, m = _full_match()
    op = create_write_operation("spec", argparse.Namespace(trim=TrimMode.BARCODES),
                                read, m, ResolutionType.FULL_MATCH)
    assert op.sequence == P1 + INSERT + P2_RC
    assert op.p1_location == (0, 3)
    assert op.p2_location == (12, 15)
    assert op.b1_location == (-4, -1)
    assert op.b2_location == (16, 19)


def test_barcode_locations_follow_selected_barcodes():
    read, m = _full_match()
    # A tied second barcode that matched one base longer; sorts after AAAA
    # only by insertion order, so a caller that resolved to it must get its
    # coordinates, not the first entry's.
    m.add_barcode_match(_hit(0, 4), "AAAAT", False, Barcode.B1)
    args = argparse.Namespace(trim=TrimMode.NONE)
    op_first = create_write_operation("specA", args, read, m, ResolutionType.DEREPLICATED_FULL,
                                      b1=B1, b2="CCCC")
    op_second = create_write_operation("specB", args, read, m, ResolutionType.DEREPLICATED_FULL,
                                       b1="AAAAT", b2="CCCC")
    assert op_first.b1_location == (0, 3)
    assert op_second.b1_location == (0, 4)
