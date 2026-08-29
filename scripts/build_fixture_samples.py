#!/usr/bin/env python3
"""Select the predeclared three-study/three-route Phase-1 smoke panel."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "config" / "samples.tsv"
OUTPUT = ROOT / "tests" / "fixtures" / "samples.tsv"
SELECTED = ["Church2025__FM-16644", "Ahuja2024__CWD1", "Ahuja2026__NA19"]

with SOURCE.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fields = reader.fieldnames
    rows = {row["sample_id"]: row for row in reader}
missing = [sample for sample in SELECTED if sample not in rows]
if missing:
    raise ValueError(f"fixture samples missing from production table: {missing}")
selected = [rows[sample] for sample in SELECTED]
if {row["study"] for row in selected} != {"Church2025", "Ahuja2024", "Ahuja2026"}:
    raise ValueError("fixture panel does not cover all three studies")
if {row["host_route"] for row in selected} != {"P_physalis", "N_septata", "none"}:
    raise ValueError("fixture panel does not cover all three host routes")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
temporary = OUTPUT.with_suffix(".tsv.tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(selected)
temporary.replace(OUTPUT)
print(f"wrote {len(selected)} fixture samples -> {OUTPUT}")
