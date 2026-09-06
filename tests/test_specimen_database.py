"""Specimens must resolve for primers that alias another primer's sequence.

Matching searches one PrimerInfo per distinct sequence and specimen lookup
compares primer objects by identity, so a specimen registered under a
second name for the same sequence must reference the canonical object.
"""

from specimux.constants import Primer
from specimux.databases import PrimerDatabase, Specimens
from specimux.models import PrimerInfo


def _db():
    db = PrimerDatabase()
    db.add_primer(PrimerInfo("fwdA", "GGGGCCCCAAAA", Primer.FWD, ["poolA"]), ["poolA"])
    db.add_primer(PrimerInfo("fwdB", "GGGGCCCCAAAA", Primer.FWD, ["poolB"]), ["poolB"])
    db.add_primer(PrimerInfo("rev", "TTTTAAAACCCC", Primer.REV, ["poolA", "poolB"]),
                  ["poolA", "poolB"])
    return db


def test_alias_specimen_resolves_through_canonical_primer():
    sp = Specimens(_db())
    sp.add_specimen("S1", "poolA", "AAAA", "fwdA", "CCCC", "rev")
    sp.add_specimen("S2", "poolB", "TTTT", "fwdB", "CCCC", "rev")

    fwd = sp.get_primers(Primer.FWD)
    assert [p.name for p in fwd] == ["fwdA"]
    rev = sp.get_primers(Primer.REV)[0]

    assert sp.specimen_for_exact_match("AAAA", "CCCC", fwd[0], rev) == "S1"
    assert sp.specimen_for_exact_match("TTTT", "CCCC", fwd[0], rev) == "S2"
    assert sp.specimens_for_barcodes_and_primers(["TTTT"], ["CCCC"], fwd[0], rev) == ["S2"]
    # The alias's pool is visible on the canonical primer
    assert set(fwd[0].pools) == {"poolA", "poolB"}
    assert set(fwd[0].barcodes) == {"AAAA", "TTTT"}


def test_alias_via_wildcard_is_deduplicated():
    sp = Specimens(_db())
    sp.add_specimen("S1", "poolA", "AAAA", "*", "CCCC", "rev")
    fwd = sp.get_primers(Primer.FWD)
    assert len(fwd) == 1
    assert sp.specimen_for_exact_match("AAAA", "CCCC", fwd[0], sp.get_primers(Primer.REV)[0]) == "S1"


def test_distinct_sequences_stay_distinct():
    db = PrimerDatabase()
    db.add_primer(PrimerInfo("f1", "GGGGCCCCAAAA", Primer.FWD, ["pool"]), ["pool"])
    db.add_primer(PrimerInfo("f2", "GGGGCCCCAAAT", Primer.FWD, ["pool"]), ["pool"])
    db.add_primer(PrimerInfo("rev", "TTTTAAAACCCC", Primer.REV, ["pool"]), ["pool"])
    sp = Specimens(db)
    sp.add_specimen("S1", "pool", "AAAA", "f1", "CCCC", "rev")
    sp.add_specimen("S2", "pool", "AAAA", "f2", "CCCC", "rev")
    f1, f2 = sp.get_primers(Primer.FWD)
    rev = sp.get_primers(Primer.REV)[0]
    assert sp.specimen_for_exact_match("AAAA", "CCCC", f1, rev) == "S1"
    assert sp.specimen_for_exact_match("AAAA", "CCCC", f2, rev) == "S2"
