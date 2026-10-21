#!/usr/bin/env python3
"""Reproduce the original Physalia/Nanomia fish input manifests."""

import argparse
import csv
import json
from pathlib import Path


def read_table(path: str) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all-host-depleted",
        action="store_true",
        help="Reproduce the frozen 161-library Physalia/Nanomia manifest; Halistemma has a separate manifest",
    )
    args = parser.parse_args()
    metadata = {
        r["analysis_sample_id"]: r
        for r in read_table("figures/cohort/sample_metadata.tsv")
    }
    positive = {
        r["analysis_sample_id"]: r["selection_basis"]
        for r in read_table("config/fish_12s_known_positive.tsv")
    }
    root = Path("data/results/phase3_cohort/inputs").resolve()
    fish, cohort = [], []
    for record in read_table("config/phase3_assembly.tsv"):
        sample = record["sample_id"]
        if args.all_host_depleted and record["host_route"] in {
            "P_physalis",
            "N_septata",
        }:
            # Preserve the original cohort-selection label in the frozen manifest.
            # Crustacean screening has moved to its own branch.
            target, basis = cohort, "both_fish_12S_and_crustacean_COI"
        elif args.all_host_depleted:
            continue
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
        if args.all_host_depleted:
            row["post_qc_pre_depletion_pairs"] = metrics["input_pairs"]
        target.append(row)
    if args.all_host_depleted:
        if (
            len(cohort) != 161
            or sum(r["host_route"] == "N_septata" for r in cohort) != 10
        ):
            raise ValueError("Unexpected host-depleted cohort membership")
        outputs = [("prey_marker_host_depleted.tsv", cohort)]
    elif len(fish) != len(positive):
        raise ValueError("Unexpected pilot cohort membership")
    else:
        outputs = [("fish_12s_host_depleted.tsv", fish)]
    for name, rows in outputs:
        with Path("config", name).open("w") as handle:
            writer = csv.DictWriter(
                handle, list(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        print(name, len(rows), "libraries;", sum(r["pairs"] for r in rows), "pairs")


if __name__ == "__main__":
    main()
