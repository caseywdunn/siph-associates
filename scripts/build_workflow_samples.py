#!/usr/bin/env python3
"""Build the normalized Snakemake sample table from the frozen manifest."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "manifest.csv"
OUTPUT = ROOT / "config" / "samples.tsv"


def safe_id(library_id: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "__", library_id)
    if not value or value in {".", ".."}:
        raise ValueError(f"unsafe library_id: {library_id!r}")
    return value


with SOURCE.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

fields = [
    "sample_id", "library_id", "specimen_id", "study", "host_route",
    "host_reference", "read_pairs", "r1_paths", "r2_paths", "include_primary",
]
output_rows = []
for row in rows:
    output_rows.append({
        "sample_id": safe_id(row["library_id"]),
        "library_id": row["library_id"],
        "specimen_id": row["specimen_id"],
        "study": row["study"],
        "host_route": row["host_reference"],
        "host_reference": row["host_reference"],
        "read_pairs": row["read_pairs"],
        "r1_paths": row["r1_paths"],
        "r2_paths": row["r2_paths"],
        "include_primary": row["include_primary"],
    })

sample_ids = [row["sample_id"] for row in output_rows]
if len(sample_ids) != len(set(sample_ids)):
    raise ValueError("safe sample IDs are not unique")
if len(output_rows) != 205:
    raise ValueError(f"expected 205 manifest rows, found {len(output_rows)}")

OUTPUT.parent.mkdir(parents=True, exist_ok=True)
temporary = OUTPUT.with_suffix(".tsv.tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(output_rows)
if OUTPUT.exists() and temporary.read_bytes() == OUTPUT.read_bytes():
    temporary.unlink()
    action = "unchanged"
else:
    temporary.replace(OUTPUT)
    action = "wrote"
print(f"{action} {len(output_rows)} workflow samples -> {OUTPUT}")
