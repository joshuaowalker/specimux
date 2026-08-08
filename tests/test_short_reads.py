"""Regression tests for reads shorter than the primer search length.

Before the fix, match_one_end computed search_start = len(sequence) - search_len,
which goes negative for short reads. align_seq treats only -1 as a sentinel, so
the negative value flowed into Python negative slicing and adjust_start(),
producing coordinates shifted by the full (search_len - len(sequence)) offset.
"""

from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from specimux.constants import Primer, Barcode
from specimux.demultiplex import match_one_end
from specimux.models import CandidateMatch, MatchParameters, PrimerInfo


PRIMER = "ACGGTCATTTAGAGGAAGTA"  # 20nt, no repeats of concern
BARCODE = "AACCGGTTAACCG"  # 13nt


def _match_short_read(read_len_padding: int) -> CandidateMatch:
    """Build a read shorter than search_len and run match_one_end on it.

    Layout (forward orientation): [barcode][primer][padding]
    match_one_end receives the reverse-complemented view, where primer_rc
    sits near the end followed by barcode_rc.
    """
    forward = BARCODE + PRIMER + "T" * read_len_padding
    rc_view = reverse_complement(forward)
    assert len(rc_view) < 80  # shorter than default search_len

    primer_info = PrimerInfo("TESTF", PRIMER, Primer.FWD, ["testpool"])
    primer_info.barcodes.add(BARCODE)
    primer_info.barcode_pairs.append((BARCODE, reverse_complement(BARCODE)))

    params = MatchParameters({PRIMER: 3}, 2, 80, False)
    seq = SeqRecord(Seq(forward), id="short_read")
    match = CandidateMatch(seq, len(BARCODE))
    match_one_end(None, match, params, rc_view, True, primer_info,
                  Primer.FWD, Barcode.B1)
    return match


def test_short_read_primer_location_is_correct():
    match = _match_short_read(read_len_padding=15)
    assert match.p1_match is not None

    location = match.get_p1_location()
    assert location is not None
    start, end = location
    # In the forward orientation of the stored read the primer occupies
    # positions len(BARCODE) .. len(BARCODE)+len(PRIMER)-1, i.e. 13..32.
    # The pre-fix code shifted coordinates by (search_len - read_len).
    assert abs(start - len(BARCODE)) <= 3
    assert abs(end - (len(BARCODE) + len(PRIMER) - 1)) <= 3


def test_short_read_barcode_location_is_correct():
    match = _match_short_read(read_len_padding=15)
    assert match.has_b1_match()
    assert BARCODE in match.best_b1()

    location = match.get_barcode1_location()
    assert location is not None
    start, end = location
    read_len = len(BARCODE) + len(PRIMER) + 15
    # Barcode occupies forward positions 0..12; pre-fix code reported
    # coordinates past the end of the read.
    assert 0 <= start <= 3
    assert end < read_len
