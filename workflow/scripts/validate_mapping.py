#!/usr/bin/env python3
"""Check Phase-5 coverage tables: complete sample x target matrices, unique keys, and accounting."""
import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--samples", required=True, nargs="+")
parser.add_argument("--bacterial", required=True)
parser.add_argument("--viral", required=True)
parser.add_argument("--libraries", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--viral-fasta", required=True)
parser.add_argument("--library-qc", required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


errors = []
samples = set(args.samples)
bacterial_targets = {r["catalog_id"] for r in rows(args.manifest) if r["catalog_id"]}
viral_targets = {line[1:].split()[0] for line in open(args.viral_fasta) if line.startswith(">")}
for label, path, targets in (("bacterial", args.bacterial, bacterial_targets), ("viral", args.viral, viral_targets)):
    table = rows(path)
    keys = Counter((r["sample_id"], r["target_id"]) for r in table)
    if any(count > 1 for count in keys.values()):
        errors.append(f"{label}: duplicated sample/target keys")
    if set(keys) != {(s, t) for s in samples for t in targets}:
        errors.append(f"{label}: matrix is not the complete sample x target product "
                      f"({len(keys)} of {len(samples) * len(targets)})")
    for r in table:
        if not 0 <= float(r["covered_fraction"]) <= 1:
            errors.append(f"{label}: covered_fraction out of range for {r['sample_id']} {r['target_id']}")
            break

qc = {r["sample_id"]: r for r in rows(args.library_qc)}
libraries = {r["sample_id"]: r for r in rows(args.libraries)}
if set(libraries) != samples:
    errors.append("library table does not match the requested samples")
counts = Counter()
for r in rows(args.bacterial):
    counts[(r["sample_id"], "bacterial")] += int(r["read_count"])
for r in rows(args.viral):
    counts[(r["sample_id"], "viral")] += int(r["read_count"])
for sample, row in libraries.items():
    pairs = int(row["input_pairs"])
    if row["read_source"] == "trimmed_regenerated":
        # Regenerated reads must reproduce the screened Phase-2 trimming exactly.
        if pairs != int(qc[sample]["pairs_after_fastp"]):
            errors.append(f"regenerated trimmed pairs differ from Phase 2: {sample}")
    for kind in ("bacterial", "viral"):
        filtered = int(row[f"{kind}_filtered_reads"])
        if filtered > 2 * pairs:
            errors.append(f"{kind} filtered reads exceed input reads: {sample}")
        if counts[(sample, kind)] > filtered:
            errors.append(f"{kind} CoverM read counts exceed filtered alignments: {sample}")

lines = [f"status\t{'PASS' if not errors else 'FAIL'}", f"samples\t{len(samples)}",
         f"bacterial_targets\t{len(bacterial_targets)}", f"viral_targets\t{len(viral_targets)}",
         f"regenerated_libraries\t{sum(r['read_source'] == 'trimmed_regenerated' for r in libraries.values())}",
         f"errors\t{len(errors)}", *[f"ERROR\t{e}" for e in errors]]
print("\n".join(lines))
if errors:
    raise SystemExit(1)
temporary = args.output.with_name(args.output.name + ".tmp")
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, args.output)
