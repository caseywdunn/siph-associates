#!/usr/bin/env python3
"""Freeze Phase-3 assembly membership from an accepted input-eligibility table."""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
FIELDS = ["sample_id", "host_route", "strategy", "retained_pairs",
          "minimum_assembly_input_pairs", "assembly_eligible"]

parser = argparse.ArgumentParser()
parser.add_argument("--fixture", action="store_true",
                    help="freeze the three-library fixture at its declared fixture floor")
args = parser.parse_args()

settings = json.loads((ROOT / "config" / "phase3_assembly.json").read_text())
if args.fixture:
    work = ROOT / "tests" / "work"
    floor = int(settings["fixture_minimum_assembly_input_pairs"])
    output = ROOT / settings["fixture_membership"]
else:
    work = ROOT / "data" / "results"
    floor = int(json.loads((ROOT / "config" / "phase3_cohort.json").read_text())["minimum_assembly_input_pairs"])
    output = ROOT / settings["membership"]

sentinel = work / "stages" / "phase3_inputs.done"
values = dict(line.split("\t", 1) for line in sentinel.read_text().splitlines() if "\t" in line)
if values.get("status") != "PASS":
    raise SystemExit(f"{sentinel} is not PASS")

with (work / "phase3_cohort" / "input_eligibility.tsv").open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len({row["sample_id"] for row in rows}) != len(rows):
    raise SystemExit("eligibility table has duplicate sample IDs")
frozen = []
for row in sorted(rows, key=lambda row: row["sample_id"]):
    eligible = int(row["retained_pairs"]) >= floor
    if not args.fixture and (row["assembly_eligible"] == "true") != eligible:
        raise SystemExit(f"eligibility disagrees with the locked floor: {row['sample_id']}")
    frozen.append({**{field: row[field] for field in FIELDS[:4]},
                   "minimum_assembly_input_pairs": floor,
                   "assembly_eligible": str(eligible).lower()})

output.parent.mkdir(parents=True, exist_ok=True)
temporary = output.with_name(output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(frozen)
os.replace(temporary, output)
print(f"{output.relative_to(ROOT)}: {len(frozen)} libraries, "
      f"{sum(row['assembly_eligible'] == 'true' for row in frozen)} eligible at {floor} pairs")
