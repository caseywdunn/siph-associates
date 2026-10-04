#!/usr/bin/env python3
"""Check eukaryote evidence grades, nonredundant reporting, and verification."""
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


settings = json.loads(Path(args.config).read_text())
gate, roles = settings["presence"], settings["role_field"]
grades = rows(args.grades)
errors = []
keys = Counter((r["sample_id"], r["reporting_unit"]) for r in grades)
if any(n > 1 for n in keys.values()):
    errors.append("duplicated sample x unit rows")
validated_units = {
    (r["sample_id"], r["reporting_unit"])
    for r in grades if r["grade"] in ("validated", "high_confidence")
}
for r in grades:
    pairs, assembled = int(r["read_pairs_ge97"]), int(r["assembled_sequences"])
    reads_ok = pairs >= gate["validated_reads"]["min_read_pairs"]
    expected = ("high_confidence" if reads_ok and assembled else "validated" if reads_ok or assembled
                else "trace")
    if r["grade"] != expected:
        errors.append(f"grade disagrees with the rule: {r['sample_id']} {r['reporting_unit']}")
    if r["kind"] in ("host", "cnidarian", "negative_control") and not assembled:
        errors.append(f"read-only detection of an excluded class: {r['sample_id']} {r['reporting_unit']}")
    allowed = {"parasite", "prey", "unassigned"} | {o["role"] for o in roles.get("verified_overrides", [])}
    if r["role"] not in allowed:
        errors.append(f"invalid role: {r['sample_id']}")
    descendants = sorted(
        unit for sample, unit in validated_units
        if sample == r["sample_id"] and unit.startswith(r["reporting_unit"] + ";")
    )
    should_count = r["grade"] in ("validated", "high_confidence") and not descendants
    if r.get("count_as_detection") != str(should_count).lower():
        errors.append(f"incorrect detection count flag: {r['sample_id']} {r['reporting_unit']}")
    expected_status = "unresolved_ancestor" if descendants else "detection" if should_count else "trace_evidence"
    if r.get("reporting_status") != expected_status:
        errors.append(f"incorrect reporting status: {r['sample_id']} {r['reporting_unit']}")
    if r.get("overlapping_descendant_units") != "|".join(descendants):
        errors.append(f"incorrect descendant annotation: {r['sample_id']} {r['reporting_unit']}")
    identity = float(r["best_identity"])
    ceiling = "genus" if identity >= 97 else "family_or_order" if identity >= 90 else "class_or_phylum"
    if r.get("naming_ceiling") != ceiling or "naming_depth" in r:
        errors.append(f"incorrect naming ceiling: {r['sample_id']} {r['reporting_unit']}")
    lca = r.get("lca_lineage", "")
    if not lca or r.get("lca_terminal_taxon") != lca.split(";")[-1]:
        errors.append(f"missing or inconsistent LCA name: {r['sample_id']} {r['reporting_unit']}")
verified = {r["sequence_id"] for r in rows(args.verification)}
missing = [r["sequence_id"] for r in rows(args.assembled) if r["sequence_id"] not in verified]
if missing:
    errors.append(f"{len(missing)} assembled non-host sequences lack NCBI verification")

evidence = [r for r in grades if r["grade"] != "trace"]
present = [r for r in grades if r.get("count_as_detection") == "true"]
lines = [f"status\t{'PASS' if not errors else 'FAIL'}",
         f"validated_evidence_rows\t{len(evidence)}",
         f"validated_detections\t{len(present)}", f"libraries\t{len({r['sample_id'] for r in present})}",
         f"high_confidence\t{sum(r['grade'] == 'high_confidence' for r in present)}",
         f"validated_ancestor_evidence_rows\t{sum(r.get('reporting_status') == 'unresolved_ancestor' for r in evidence)}",
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
