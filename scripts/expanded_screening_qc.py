#!/usr/bin/env python3
"""Extend the original screen inventory using completed Bargmannia products."""

import csv
import json
from pathlib import Path


def read(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    config = json.loads(Path("config/reporting_sources.json").read_text())
    rows = read("data/results/aggregation/cohort/library_qc.tsv")
    fields = ["sample_id", "status", "bracken_genus_rows", "bracken_species_rows"]
    root = Path("data/results/bargmannia/2026-10-19")
    if "status\tPASS" not in (root / "stages/phase3_assembly.done").read_text():
        raise ValueError("Bargmannia workflow must have passed before reporting")
    for sample in read("data/metadata/manko2026/sample_metadata.tsv"):
        sid = sample["analysis_sample_id"]
        row = {field: sample.get(field, "") for field in fields}
        row.update(
            sample_id=sid,
            raw_pairs=sample["read_pairs"],
            status="PASS",
            bracken_genus_rows=len(read(root / "screens/bracken" / f"{sid}.G.tsv")),
            bracken_species_rows=len(read(root / "screens/bracken" / f"{sid}.S.tsv")),
        )
        rows.append(row)
    expected = {r["sample_id"] for r in read(config["samples"])}
    if len(rows) != len(expected) or {r["sample_id"] for r in rows} != expected:
        raise ValueError("Screen inventory differs from the reporting cohort")
    output = Path(config["screen_qc"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
