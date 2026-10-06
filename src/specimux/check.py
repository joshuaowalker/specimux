#!/usr/bin/env python3

"""
Validation-only mode (specimux --check).

Loads a primers file and a specimen file exactly as a demultiplexing run does,
but collects every problem instead of stopping at the first, and needs no
sequence file. Intended for rejecting bad inputs before sequencing data is
uploaded or basecalled.
"""

import json
from dataclasses import asdict, dataclass, field
from typing import List

from .io_utils import load_primers, load_specimens
from .models import InputProblem


@dataclass
class CheckResult:
    """Outcome of checking a primers file and a specimen file."""
    problems: List[InputProblem] = field(default_factory=list)
    primers: int = 0
    pools: int = 0
    specimens: int = 0

    @property
    def valid(self) -> bool:
        return not self.problems


def check_inputs(primer_file: str, specimen_file: str) -> CheckResult:
    """Run every input check a demultiplexing run makes before reading sequences."""
    result = CheckResult()
    registry = load_primers(primer_file, result.problems)
    specimens = load_specimens(specimen_file, registry, result.problems)
    if result.valid:
        # Warns about suspicious barcodes (logged, not failures), as a real run does
        specimens.validate()

    stats = registry.get_pool_stats()
    result.primers = stats['total_primers']
    result.pools = stats['total_pools']
    result.specimens = len(specimens._specimens)
    return result


def format_text(result: CheckResult) -> str:
    """One line per problem, then a summary line."""
    lines = [str(p) for p in result.problems]
    if result.valid:
        lines.append(f"OK: {result.primers} primers in {result.pools} pools, "
                     f"{result.specimens} specimens")
    else:
        n = len(result.problems)
        lines.append(f"FAILED: {n} problem{'s' if n != 1 else ''} found")
    return "\n".join(lines)


def format_json(result: CheckResult) -> str:
    """Machine-readable result for orchestration tools."""
    return json.dumps({
        "valid": result.valid,
        "problems": [asdict(p) for p in result.problems],
        "primers": result.primers,
        "pools": result.pools,
        "specimens": result.specimens,
    }, indent=2)
