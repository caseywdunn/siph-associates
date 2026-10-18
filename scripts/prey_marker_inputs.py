#!/usr/bin/env python3
"""Freeze existing host-depleted pairs as inputs to the revised prey pilots."""

import csv
import json
from pathlib import Path


def read_table(path: str) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    metadata = {
        r["analysis_sample_id"]: r
        for r in read_table("figures/cohort/sample_metadata.tsv")
    }
    positive = {
        r["analysis_sample_id"]: r["selection_basis"]
        for r in read_table("config/fish_12s_known_positive.tsv")
    }
    root = Path("data/results/phase3_cohort/inputs").resolve()
    fish, crustacean = [], []
    for record in read_table("config/phase3_assembly.tsv"):
        sample = record["sample_id"]
        if record["host_route"] == "N_septata":
            target, basis = crustacean, "all_N_septata_libraries"
        elif record["host_route"] == "P_physalis" and sample in positive:
            target, basis = fish, positive[sample]
        else:
            continue
        metrics = json.loads((root / f"{sample}.json").read_text())
        assert metrics["strategy"] == record["strategy"] == "reference_depleted"
        assert metrics["retained_pairs"] == int(record["retained_pairs"])
        row = {
            "analysis_sample_id": sample,
            "host": metadata[sample]["species_current"],
            "pairs": metrics["retained_pairs"],
            "read_source": "reference_depleted",
            "selection_basis": basis + ";all_retained_pairs;no_additional_subsampling",
            "read_pair_id": "host_depleted",
            "host_route": record["host_route"],
            "upstream_metrics": str(root / f"{sample}.json"),
        }
        for mate in (1, 2):
            path = root / f"{sample}_R{mate}.fastq.gz"
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f"Missing retained input: {path}")
            row[f"r{mate}_path"] = str(path)
            row[f"r{mate}_bytes"] = path.stat().st_size
        target.append(row)
    if len(fish) != len(positive) or len(crustacean) != 10:
        raise ValueError("Unexpected pilot cohort membership")
    for name, rows in [
        ("fish_12s_host_depleted.tsv", fish),
        ("crustacean_coi_pilot.tsv", crustacean),
    ]:
        with Path("config", name).open("w") as handle:
            writer = csv.DictWriter(
                handle, list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        print(name, len(rows), "libraries;", sum(r["pairs"] for r in rows), "pairs")


if __name__ == "__main__":
    main()
