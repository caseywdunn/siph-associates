#!/usr/bin/env python3
"""Check Phase-6 outputs: complete grade matrices, valid grades, and non-empty analysis tables."""
import argparse
import csv
import os
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--bacterial", required=True)
parser.add_argument("--viral", required=True)
parser.add_argument("--contamination", required=True)
parser.add_argument("--primary", required=True, nargs="+")
parser.add_argument("--sensitivity", required=True, nargs="+")
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


errors = []
for label, path, allowed in (("bacterial", args.bacterial, {"none", "trace", "validated", "high_confidence"}),
                             ("viral", args.viral, {"none", "partial", "validated", "high_confidence"})):
    table = rows(path)
    keys = Counter((r["sample_id"], r["target_id"]) for r in table)
    samples, targets = {s for s, _ in keys}, {t for _, t in keys}
    if any(n > 1 for n in keys.values()) or len(keys) != len(samples) * len(targets):
        errors.append(f"{label}: grades are not a complete, unique sample x target matrix")
    if len(samples) != 205:
        errors.append(f"{label}: expected 205 libraries, found {len(samples)}")
    bad = {r["grade"] for r in table} - allowed
    if bad:
        errors.append(f"{label}: unexpected grades {sorted(bad)}")
    if label == "bacterial":
        decoy_calls = [r for r in table if r["role"] == "decoy" and r["grade"] in ("validated", "high_confidence")]
        if decoy_calls:
            errors.append(f"bacterial: {len(decoy_calls)} decoy presences at the primary threshold")
for r in rows(args.contamination):
    if r["p_value"] and not 0 < float(r["p_value"]) <= 1:
        errors.append(f"contamination p-value out of range: {r['target_id']}")
for path in args.primary + args.sensitivity:
    if not rows(path) and "phage_host_links" not in path:
        errors.append(f"empty analysis table: {Path(path).name}")

lines = [f"status\t{'PASS' if not errors else 'FAIL'}", f"errors\t{len(errors)}", *[f"ERROR\t{e}" for e in errors]]
print("\n".join(lines))
if errors:
    raise SystemExit(1)
temporary = args.output.with_name(args.output.name + ".tmp")
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, args.output)
