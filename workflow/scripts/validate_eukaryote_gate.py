#!/usr/bin/env python3
"""Check eukaryote grades against the locked rules and that every assembled sequence was verified."""
import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--grades", required=True)
parser.add_argument("--assembled", required=True)
parser.add_argument("--verification", required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


gate = json.loads(Path(args.config).read_text())["presence"]
grades = rows(args.grades)
errors = []
keys = Counter((r["sample_id"], r["reporting_unit"]) for r in grades)
if any(n > 1 for n in keys.values()):
    errors.append("duplicated sample x unit rows")
for r in grades:
    pairs, assembled = int(r["read_pairs_ge97"]), int(r["assembled_sequences"])
    reads_ok = pairs >= gate["validated_reads"]["min_read_pairs"]
    expected = ("high_confidence" if reads_ok and assembled else "validated" if reads_ok or assembled
                else "trace")
    if r["grade"] != expected:
        errors.append(f"grade disagrees with the rule: {r['sample_id']} {r['reporting_unit']}")
    if r["kind"] in ("host", "cnidarian", "negative_control") and not assembled:
        errors.append(f"read-only detection of an excluded class: {r['sample_id']} {r['reporting_unit']}")
    if r["role"] not in ("parasite", "prey", "unassigned"):
        errors.append(f"invalid role: {r['sample_id']}")
verified = {r["sequence_id"] for r in rows(args.verification)}
missing = [r["sequence_id"] for r in rows(args.assembled) if r["sequence_id"] not in verified]
if missing:
    errors.append(f"{len(missing)} assembled non-host sequences lack NCBI verification")

present = [r for r in grades if r["grade"] != "trace"]
lines = [f"status\t{'PASS' if not errors else 'FAIL'}",
         f"validated_detections\t{len(present)}", f"libraries\t{len({r['sample_id'] for r in present})}",
         f"high_confidence\t{sum(r['grade'] == 'high_confidence' for r in grades)}",
         f"trace\t{sum(r['grade'] == 'trace' for r in grades)}",
         *[f"role_{k}\t{v}" for k, v in sorted(Counter(r["role"] for r in present).items())],
         f"contamination_flags\t{sum(r['contamination'] == 'true' for r in grades)}",
         f"errors\t{len(errors)}", *[f"ERROR\t{e}" for e in errors]]
print("\n".join(lines))
if errors:
    raise SystemExit(1)
temporary = args.output.with_name(args.output.name + ".tmp")
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, args.output)
