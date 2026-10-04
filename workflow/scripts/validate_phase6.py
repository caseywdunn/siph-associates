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
parser.add_argument("--metadata", required=True)
parser.add_argument("--nanomia", required=True)
parser.add_argument("--thresholds", required=True, nargs="+")
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
for path in args.primary + args.sensitivity + args.thresholds:
    if not rows(path) and not any(name in path for name in ("phage_host_links", "nanomia_contamination", "physalia_models", "physalia_region_permutation")):
        errors.append(f"empty analysis table: {Path(path).name}")

# Explicit selection and full-threshold checks protect the corrected batch tests.
metadata = rows(args.metadata)
if len(metadata) != 205 or len({row["sample_id"] for row in metadata}) != 205:
    errors.append("flowcell metadata must retain all 205 unique libraries")
for row in metadata:
    flowcells = row["flowcells"].split(";") if row["flowcells"] else []
    eligible = row["inference_eligible"] == "true"
    if (eligible != (len(flowcells) == 1) or int(row["flowcell_count"]) != len(flowcells)
            or any(len(flowcell.split(":")) != 3 for flowcell in flowcells)
            or (eligible and row["flowcell"] != flowcells[0])):
        errors.append(f"invalid physical flowcell normalization: {row['sample_id']}")
meta = {row["sample_id"]: row for row in metadata}
bacterial = rows(args.bacterial)
eligible_physalia = {sample for sample, row in meta.items()
                    if row["inference_eligible"] == "true"
                    and row["host_species"].startswith("Physalia") and row["ocean_region"]}
primary_dir = Path(args.nanomia).parent
for threshold, directory, column in (
    (10, primary_dir, "grade"),
    (5, primary_dir.parent / "sensitivity" / "breadth_5pct", "grade_at_5pct"),
    (20, primary_dir.parent / "sensitivity" / "breadth_20pct", "grade_at_20pct"),
):
    calls = [row for row in bacterial if row["role"] != "decoy"
             and row[column] in {"validated", "high_confidence"}]
    expected = Counter(row["target_id"] for row in calls)
    incidence = rows(directory / "bacterial_incidence.tsv")
    for row in incidence:
        if int(row["present_libraries"]) != expected[row["target_id"]]:
            errors.append(f"{threshold}% incidence disagrees with grades: {row['target_id']}")
    designs = rows(directory / "physalia_design.tsv")
    if sum(int(row["N"]) for row in designs) != len(eligible_physalia):
        errors.append(f"{threshold}% Physalia design has an incorrect eligible sample count")
    physalia_counts = Counter(row["target_id"] for row in calls if row["sample_id"] in eligible_physalia)
    for row in rows(directory / "physalia_region_permutation.tsv"):
        if (int(row["presences"]) != physalia_counts[row["target_id"]]
                or int(row["libraries"]) != len(eligible_physalia)
                or row["grade_column"] != column):
            errors.append(f"{threshold}% Physalia inference uses the wrong sample set or threshold")
    for name in ("contamination_tests", "nanomia_contamination_tests"):
        path = Path(args.contamination) if threshold == 10 and name == "contamination_tests" else directory / f"{name}.tsv"
        for row in rows(path):
            target_calls = [call for call in calls if call["target_id"] == row["target_id"]
                            and (name == "contamination_tests" or
                                 meta[call["sample_id"]]["host_species"].startswith("Nanomia"))]
            n_eligible = sum(meta[call["sample_id"]]["inference_eligible"] == "true" for call in target_calls)
            if (int(row["validated_presences"]) != len(target_calls)
                    or int(row["eligible_presences"]) != n_eligible or row["grade_column"] != column):
                errors.append(f"{threshold}% {name} uses the wrong sample set or threshold")
            if row["status"] == "tested" and not (0 < float(row["p_value"]) <= 1
                                                   and 0 < float(row["q_value"]) <= 1):
                errors.append(f"{threshold}% {name}: invalid p/q value")

lines = [f"status\t{'PASS' if not errors else 'FAIL'}", f"errors\t{len(errors)}", *[f"ERROR\t{e}" for e in errors]]
print("\n".join(lines))
if errors:
    raise SystemExit(1)
temporary = args.output.with_name(args.output.name + ".tmp")
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, args.output)
