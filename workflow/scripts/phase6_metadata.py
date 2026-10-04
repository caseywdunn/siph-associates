#!/usr/bin/env python3
"""Normalize physical flowcells and make the inferential sample selection explicit."""

import argparse
import csv
from pathlib import Path


def physical_flowcells(batches: str) -> tuple[str, ...]:
    """Collapse lanes using manifest instrument:run:flowcell:lane identifiers."""
    flowcells = set()
    for batch in batches.split(";"):
        if not batch:
            continue
        fields = batch.split(":")
        if len(fields) != 4 or not all(fields):
            raise ValueError(f"Expected instrument:run:flowcell:lane, got {batch!r}")
        flowcells.add(":".join(fields[:3]))
    return tuple(sorted(flowcells))


def normalize_manifest(path: Path) -> list[dict[str, str | int]]:
    """Retain every primary library; exclude ambiguous batches only from inference."""
    records = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["include_primary"].lower() != "true":
                continue
            flowcells = physical_flowcells(row["sequencing_batches"])
            reason = "" if len(flowcells) == 1 else (
                "multiple_flowcells" if flowcells else "missing_flowcell"
            )
            records.append({
                "sample_id": f"{row['study']}__{row['library_id'].split(':')[-1]}",
                "study": row["study"],
                "host_species": row["species_current"],
                "ocean_region": row["ocean_region"],
                "flowcells": ";".join(flowcells),
                "flowcell_count": len(flowcells),
                "inference_policy": "single_flowcell",
                "flowcell": flowcells[0] if len(flowcells) == 1 else "",
                "inference_eligible": str(len(flowcells) == 1).lower(),
                "inference_exclusion": reason,
            })
    return sorted(records, key=lambda row: row["sample_id"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--flowcell-policy", required=True, choices=("single_flowcell",))
    args = parser.parse_args()
    records = normalize_manifest(args.manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    temporary.replace(args.output)
    print(f"libraries={len(records)} single_flowcell="
          f"{sum(row['inference_eligible'] == 'true' for row in records)}")


if __name__ == "__main__":
    main()
