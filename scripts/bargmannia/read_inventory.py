#!/usr/bin/env python3
"""Validate fastp accounting and provide the shared fish screen's input table."""

import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.samples.open() as handle:
        rows = [
            r
            for r in csv.DictReader(handle, delimiter="\t")
            if r["analysis_sample_id"] == args.sample
        ]
    assert len(rows) == 1
    sample = rows[0]
    provenance = json.loads(
        (args.root / "provenance/trim" / f"{args.sample}.json").read_text()
    )
    fastp = json.loads((args.root / "trim" / f"{args.sample}.fastp.json").read_text())
    assert fastp["summary"]["before_filtering"]["total_reads"] == 2 * min(
        int(sample["read_pairs"]), provenance["parameters"]["cap_pairs"]
    )
    assert (
        fastp["summary"]["after_filtering"]["total_reads"]
        == 2 * provenance["pairs_after_fastp"]
    )
    r1 = args.root / "retained/trimmed" / f"{args.sample}_R1.fastq.gz"
    r2 = args.root / "retained/trimmed" / f"{args.sample}_R2.fastq.gz"
    row = {
        "analysis_sample_id": args.sample,
        "host": sample["species_current"],
        "read_pair_id": "capped_trimmed",
        "r1_path": str(r1),
        "r2_path": str(r2),
        "r1_bytes": r1.stat().st_size,
        "r2_bytes": r2.stat().st_size,
        "pairs": provenance["pairs_after_fastp"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
