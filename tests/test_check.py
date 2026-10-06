"""Validation-only mode: specimux --check primers.fasta specimens.txt.

The check must reject everything a real run rejects from these two files,
report every problem at once, and need no sequence file.
"""

import json
import subprocess
import sys

import pytest

from specimux.check import check_inputs
from specimux.io_utils import read_primers_file, read_specimen_file

HEADER = "SampleID\tPrimerPool\tFwIndex\tFwPrimer\tRvIndex\tRvPrimer\n"


def _write(path, text):
    path.write_text(text)
    return str(path)


def _run_check(*args):
    return subprocess.run([sys.executable, "-c", "from specimux.cli import main; main()",
                           "--check", *args], capture_output=True, text=True)


@pytest.fixture
def run155(temp_dir):
    """Primers and index shaped like the Run155 failure: pool ITS has only a
    reverse primer, and every specimen names a forward primer that does not exist."""
    primers = _write(temp_dir / "primers.fasta",
                     ">ITS1Fngs pool=ITSnew position=forward\nGGTCATTTAGAGGAAGTAA\n"
                     ">ITS4ngsUni pool=ITS,ITSnew position=reverse\nCCTSCSCTTANTDATATGC\n")
    rows = "".join(f"S{i}\tITS\tAAAA{i:04d}\tITS1F\tCCCC{i:04d}\tITS4ngsUni\n" for i in range(50))
    specimens = _write(temp_dir / "Index.txt", HEADER + rows)
    return primers, specimens


def test_run155_reports_both_problems(run155):
    primers, specimens = run155
    result = check_inputs(primers, specimens)
    assert not result.valid
    assert [(p.file, p.line, p.count) for p in result.problems] == [
        (primers, 3, 1), (specimens, 2, 50)]
    assert "Pool ITS has no forward primers" in result.problems[0].message
    assert "'ITS1F' is not in the primers file" in result.problems[1].message
    assert "ITS1Fngs" in result.problems[1].message


def test_valid_inputs(integration_test_data):
    result = check_inputs(str(integration_test_data / "primers.fasta"),
                          str(integration_test_data / "specimens.txt"))
    assert result.valid
    assert result.specimens > 0


def test_every_problem_is_reported(temp_dir):
    primers = _write(temp_dir / "p.fasta",
                     ">A pool=P1 position=forward\nACGTACGTACGTAC\n"
                     ">B pool=P1 position=sideways\nACGTACGTACGTAA\n"
                     ">C position=reverse\nACGTACGTACGTAG\n"
                     ">A pool=P1 position=forward\nACGTACGTACGTCC\n"
                     ">R pool=P1 position=reverse\nTTTTACGTACGTCC\n")
    specimens = _write(temp_dir / "s.tsv", HEADER +
                       "S1\tP1\tAAAA\tA\tCCCC\tR\n"
                       "S1\tP1\tAAAT\tA\tCCCA\tR\n"
                       "S2\tP9\tAAAG\tA\tCCCG\tR\n"
                       "S3\tP1\t\tR\tCCCT\tA\n"
                       "S4\tP1\tAAGG\n")
    found = [(p.line, p.message) for p in check_inputs(primers, specimens).problems]
    assert found == [
        (3, "Invalid primer position 'sideways' for B (expected forward or reverse)"),
        (5, "Missing pool specification for primer C"),
        (7, "Duplicate primer name: A"),
        (3, "Duplicate SampleID 'S1'"),
        (4, "PrimerPool 'P9' is not defined in the primers file"),
        (5, "FwIndex is empty (single-indexed demultiplexing is not supported)"),
        (5, "FwPrimer 'R' is not a forward primer in the primers file"),
        (5, "RvPrimer 'A' is not a reverse primer in the primers file"),
        (6, "Row has 3 fields but the header has 6"),
    ]


def test_missing_columns_and_unreadable_files(temp_dir):
    specimens = _write(temp_dir / "s.csv", HEADER.replace("\t", ","))
    problems = check_inputs(str(temp_dir / "missing.fasta"), specimens).problems
    assert "Cannot read primers file" in problems[0].message
    assert problems[1].line == 1 and "Missing required columns" in problems[1].message


def test_real_run_loaders_reject_the_same_inputs(run155):
    primers, specimens = run155
    with pytest.raises(ValueError, match="Pool ITS has no forward primers"):
        read_primers_file(primers)


def test_real_run_specimen_loader_lists_all_problems(temp_dir):
    primers = _write(temp_dir / "p.fasta",
                     ">F pool=P position=forward\nACGTACGTACGTAC\n"
                     ">R pool=P position=reverse\nTTTTACGTACGTCC\n")
    specimens = _write(temp_dir / "s.tsv", HEADER +
                       "S1\tP\tAAAA\tX\tCCCC\tR\n"
                       "S2\tP\tAAAT\tF\tCCCA\tY\n")
    with pytest.raises(ValueError) as e:
        read_specimen_file(specimens, read_primers_file(primers))
    assert "FwPrimer 'X'" in str(e.value) and "RvPrimer 'Y'" in str(e.value)


def test_cli_exit_codes_and_json(run155, integration_test_data):
    ok = _run_check(str(integration_test_data / "primers.fasta"),
                    str(integration_test_data / "specimens.txt"))
    assert ok.returncode == 0, ok.stderr
    assert ok.stdout.startswith("OK:")

    bad = _run_check(*run155)
    assert bad.returncode == 1
    assert bad.stdout.splitlines()[-1] == "FAILED: 2 problems found"

    as_json = _run_check("--json", *run155)
    assert as_json.returncode == 1
    report = json.loads(as_json.stdout)
    assert report["valid"] is False
    assert [p["line"] for p in report["problems"]] == [3, 2]


def test_sequence_file_still_required_without_check(integration_test_data):
    result = subprocess.run([sys.executable, "-c", "from specimux.cli import main; main()",
                             str(integration_test_data / "primers.fasta"),
                             str(integration_test_data / "specimens.txt")],
                            capture_output=True, text=True)
    assert result.returncode == 2
    assert "sequence_file" in result.stderr
