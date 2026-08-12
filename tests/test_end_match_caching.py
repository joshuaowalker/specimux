"""Tests for per-(primer, strand) end-match caching.

The cache shares EndMatch results across CandidateMatches, so assembly
must hand each match independent AlignmentResult copies: trim_locations()
mutates stored locations in place, and a shared object would corrupt the
coordinates of every other match built from the same cached end.
"""

from Bio.Seq import reverse_complement
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from specimux.constants import Primer, Barcode
from specimux.demultiplex import apply_end_match, compute_end_match
from specimux.models import CandidateMatch, MatchParameters, PrimerInfo


PRIMER = "ACGGTCATTTAGAGGAAGTA"
BARCODE = "AACCGGTTAACCG"


def _end_match():
    forward = BARCODE + PRIMER + "T" * 40
    rc_view = reverse_complement(forward)
    primer_info = PrimerInfo("TESTF", PRIMER, Primer.FWD, ["testpool"])
    primer_info.barcodes.add(BARCODE)
    primer_info.barcode_pairs.append((BARCODE, reverse_complement(BARCODE)))
    params = MatchParameters({PRIMER: 3}, 2, 80, False)
    end = compute_end_match(None, params, rc_view, primer_info, Primer.FWD, Barcode.B1)
    seq = SeqRecord(Seq(forward), id="read")
    return end, primer_info, seq


def test_end_match_found():
    end, _, _ = _end_match()
    assert end.primer_match is not None
    assert len(end.barcode_matches) == 1


def test_shared_end_match_is_isolated_between_candidate_matches():
    end, primer_info, seq = _end_match()

    m1 = CandidateMatch(seq, len(BARCODE))
    m2 = CandidateMatch(seq, len(BARCODE))
    apply_end_match(m1, end, primer_info, True, Primer.FWD, Barcode.B1)
    apply_end_match(m2, end, primer_info, True, Primer.FWD, Barcode.B1)

    p1_before = m2.get_p1_location()
    b1_before = m2.get_barcode1_location()

    # Mutating one match's coordinates must not leak into the other,
    # nor into the cached EndMatch itself.
    m1.trim_locations(10)

    assert m2.get_p1_location() == p1_before
    assert m2.get_barcode1_location() == b1_before

    m3 = CandidateMatch(seq, len(BARCODE))
    apply_end_match(m3, end, primer_info, True, Primer.FWD, Barcode.B1)
    assert m3.get_p1_location() == p1_before
    assert m3.get_barcode1_location() == b1_before
