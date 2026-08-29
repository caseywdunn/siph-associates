#!/usr/bin/env python3
"""Derive the Phase-0 storage/scale estimate from frozen metadata."""
from __future__ import annotations

import csv
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CAP = 400_000_000
manifest = list(csv.DictReader((REPO / "manifest.csv").open()))
raw = list(csv.DictReader(
    (REPO / "data" / "metadata" / "raw_files.tsv").open(), delimiter="\t"
))
bytes_by_library = {row["library_id"]: 0 for row in manifest}
for row in raw:
    bytes_by_library[row["library_id"]] += int(row["r1_bytes"]) + int(row["r2_bytes"])

raw_pairs = sum(int(row["read_pairs"]) for row in manifest)
capped_pairs = sum(min(int(row["read_pairs"]), CAP) for row in manifest)
raw_bytes = sum(bytes_by_library.values())
capped_bytes = sum(
    bytes_by_library[row["library_id"]] *
    min(int(row["read_pairs"]), CAP) / int(row["read_pairs"])
    for row in manifest
)
capped_libraries = sum(int(row["read_pairs"]) > CAP for row in manifest)
rows = [
    ("analytical_libraries", len(manifest), "count", "deduplicated manifest"),
    ("paired_fastq_files", len(raw), "pairs", "selected original files"),
    ("raw_read_pairs", raw_pairs, "read_pairs", "sum of exact manifest depths"),
    ("raw_compressed_bytes", raw_bytes, "bytes", "sum of selected R1+R2 sizes"),
    ("libraries_above_400m_cap", capped_libraries, "count", "deterministic cap candidates"),
    ("capped_read_pairs", capped_pairs, "read_pairs", "sum min(depth, 400M)"),
    ("estimated_capped_compressed_bytes", round(capped_bytes), "bytes",
     "linear estimate from each library's retained pair fraction"),
    ("scratch_reservation_with_20pct_margin", round(capped_bytes * 1.2), "bytes",
     "trimmed-read planning estimate; actual fastp compression may differ"),
]
out = REPO / "data" / "metadata" / "storage_estimate.tsv"
with out.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
    writer.writerow(["metric", "value", "unit", "basis"])
    writer.writerows(rows)
print(f"wrote {len(rows)} metrics -> {out}")
