#!/usr/bin/env python3
"""Check that the assembly summary covers exactly the frozen membership."""
import argparse
import csv
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--summary", required=True)
parser.add_argument("--membership", required=True)
parser.add_argument("--eligibility", required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


rows = table(args.summary)
membership = {row["sample_id"]: row for row in table(args.membership)}
eligibility = {row["sample_id"]: row for row in table(args.eligibility)}
errors = []
if len(rows) != len(membership) or {row["sample_id"] for row in rows} != set(membership):
    errors.append("summary samples do not exactly match the frozen membership")
if set(membership) != set(eligibility):
    errors.append("frozen membership samples do not match the input-eligibility table")
for sample, member in membership.items():
    if member["retained_pairs"] != eligibility.get(sample, {}).get("retained_pairs"):
        errors.append(f"frozen retained pairs differ from the eligibility table: {sample}")
    eligible = int(member["retained_pairs"]) >= int(member["minimum_assembly_input_pairs"])
    if (member["assembly_eligible"] == "true") != eligible:
        errors.append(f"frozen eligibility disagrees with its floor: {sample}")
for row in rows:
    sample, status = row["sample_id"], row["assembly_status"]
    eligible = membership.get(sample, {}).get("assembly_eligible") == "true"
    if status not in (("assembled", "no_contigs") if eligible else ("excluded_below_floor",)):
        errors.append(f"unexpected assembly status {status!r}: {sample}")
        continue
    if not eligible:
        continue
    if row["checkm2_status"] not in ("assessed", "no_annotations", "no_bins"):
        errors.append(f"invalid CheckM2 status: {sample}")
    if (row["checkm2_status"] == "no_bins") != (int(row["raw_bins"]) == 0):
        errors.append(f"CheckM2 status disagrees with the bin count: {sample}")
    if int(row["assessed_bins"]) + int(row["unassessed_bins"]) != int(row["raw_bins"]):
        errors.append(f"bin accounting mismatch: {sample}")
    if row["checkm2_status"] == "assessed" and int(row["unassessed_bins"]):
        errors.append(f"CheckM2 omitted bins without recording them as unassessable: {sample}")
    if int(row["raw_bins"]) and not int(row["binnable_contigs"]):
        errors.append(f"bins reported without binnable contigs: {sample}")

count = lambda status: sum(row["assembly_status"] == status for row in rows)
lines = [f"status\t{'PASS' if not errors else 'FAIL'}", f"samples\t{len(rows)}",
         f"assembled\t{count('assembled')}", f"no_contigs\t{count('no_contigs')}",
         f"excluded_below_floor\t{count('excluded_below_floor')}",
         f"raw_bins\t{sum(int(row['raw_bins'] or 0) for row in rows)}",
         f"viral_contigs\t{sum(int(row['viral_contigs'] or 0) for row in rows)}",
         f"errors\t{len(errors)}", *[f"ERROR\t{error}" for error in errors]]
print("\n".join(lines))
if errors:
    raise SystemExit(1)
temporary = args.output.with_name(args.output.name + ".tmp")
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, args.output)
