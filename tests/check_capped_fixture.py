#!/usr/bin/env python3
"""Check that the forced capped fixture preserves mates and records the exact cap."""
import gzip
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SCRATCH = ROOT / "tests" / "cap_scratch" / "trimmed"
PROVENANCE = ROOT / "tests" / "cap_work" / "provenance" / "trim" / "Church2025__FM-16644.json"


def reads(path):
    with gzip.open(path, "rt") as handle:
        lines = sum(1 for _ in handle)
    if lines % 4:
        raise ValueError(f"incomplete FASTQ: {path}")
    return lines // 4


counts = [reads(SCRATCH / f"Church2025__FM-16644_R{mate}.fastq.gz") for mate in (1, 2)]
record = json.loads(PROVENANCE.read_text())
if counts[0] != counts[1] or counts[0] <= 0 or counts[0] > 1000:
    raise AssertionError(f"invalid post-fastp paired counts: {counts}")
if record["pairs_entering_fastp"] != 1000:
    raise AssertionError("capped fixture did not send exactly 1000 pairs to fastp")
if record["parameters"]["cap_method"] != "bbtools_exact_seeded_pair_sampling":
    raise AssertionError("unexpected cap method")
print(f"status\tPASS\npairs_entering_fastp\t1000\npairs_after_fastp\t{counts[0]}\n")
