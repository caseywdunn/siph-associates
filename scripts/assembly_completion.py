#!/usr/bin/env python3
"""Freeze remaining assembly inputs and validate existing host-depletion records."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def table(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path: str | Path, rows: list[dict], fields: list[str]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def unique(rows: list[dict], key: str) -> dict[str, dict]:
    result = {r[key]: r for r in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate {key}")
    return result


def plan(settings: dict) -> tuple[list[dict], list[dict]]:
    if (
        settings["minimum_assembly_input_pairs"] != 0
        or not settings["empty_assembly_is_valid"]
    ):
        raise ValueError(
            "Current policy requires no read-count exclusion and accepts empty assemblies"
        )
    original = unique(table(settings["original_samples"]), "sample_id")
    old_membership = unique(table(settings["original_membership"]), "sample_id")
    halistemma = unique(table(settings["halistemma_samples"]), "sample_id")
    new = unique(table(settings["bargmannia_samples"]), "analysis_sample_id")
    if (
        set(original) != set(old_membership)
        or not set(halistemma) <= set(original)
        or set(new) & set(original)
    ):
        raise ValueError("Original, replacement, and new library identities disagree")
    metadata = unique(
        table(settings["original_metadata"]) + table(settings["bargmannia_metadata"]),
        "analysis_sample_id",
    )
    completion, cohort = [], []
    for sample, source in sorted((original | new).items()):
        reason = "reuse_completed_assembly"
        input_root = settings["original_root"]
        result_root = input_root
        if sample in new:
            input_root = result_root = settings["bargmannia_root"]
        elif sample in halistemma:
            source = halistemma[sample]
            input_root = settings["halistemma_root"]
            result_root = settings["work_root"]
            reason = "replace_full_read_assembly_with_host_depleted_assembly"
        elif old_membership[sample]["assembly_eligible"] == "false":
            result_root = settings["work_root"]
            reason = "assemble_previously_excluded_library"
        record = json.loads(
            (Path(input_root) / "phase3_cohort/inputs" / f"{sample}.json").read_text()
        )
        if (
            record["sample_id"] != sample
            or record["host_route"] != source["host_route"]
        ):
            raise ValueError(f"Input identity mismatch for {sample}")
        row = {
            "sample_id": sample,
            "species_current": metadata[sample]["species_current"],
            "study": source["study"],
            "host_route": source["host_route"],
            "strategy": record["strategy"],
            "retained_pairs": record["retained_pairs"],
            "minimum_assembly_input_pairs": 0,
            "assembly_eligible": "true",
            "input_root": input_root,
            "result_root": result_root,
            "selection_reason": reason,
            "r1_paths": source["r1_paths"],
            "r2_paths": source["r2_paths"],
            "raw_pairs": source["read_pairs"],
        }
        cohort.append(row)
        if result_root == settings["work_root"]:
            completion.append(row)
    if (
        len(completion) != settings["expected_completion_samples"]
        or len(cohort) != settings["expected_cohort_samples"]
    ):
        raise ValueError("Unexpected completion/cohort library counts")
    return completion, cohort


def validate_inputs(rows: list[dict]) -> dict:
    records = []
    for row in rows:
        sample, root = row["sample_id"], Path(row["input_root"])
        metrics_path = root / "phase3_cohort/inputs" / f"{sample}.json"
        provenance_path = root / "provenance/phase3_cohort/inputs" / f"{sample}.json"
        trim_path = root / "provenance/trim" / f"{sample}.json"
        metrics, provenance, trim = [
            json.loads(p.read_text())
            for p in (metrics_path, provenance_path, trim_path)
        ]
        if any(r["sample_id"] != sample for r in (metrics, provenance, trim)):
            raise ValueError(f"Sample identity mismatch: {sample}")
        if (
            metrics["strategy"] != "reference_depleted"
            or provenance["strategy"] != "reference_depleted"
        ):
            raise ValueError(f"Unexpected assembly input strategy: {sample}")
        if any(r["host_route"] != row["host_route"] for r in (metrics, provenance)):
            raise ValueError(f"Host route mismatch: {sample}")
        if metrics["input_pairs"] != trim["pairs_after_fastp"] or provenance[
            "upstream_provenance_sha256"
        ] != sha256(trim_path):
            raise ValueError(f"Input provenance mismatch: {sample}")
        if (
            metrics["retained_pairs"] != int(row["retained_pairs"])
            or not 0 <= metrics["retained_pairs"] <= metrics["input_pairs"]
        ):
            raise ValueError(f"Invalid retained pairs: {sample}")
        for mate, output in zip((1, 2), provenance["outputs"], strict=True):
            path = root / "phase3_cohort/inputs" / f"{sample}_R{mate}.fastq.gz"
            if (
                path.resolve() != Path(output["path"]).resolve()
                or path.stat().st_size != output["size_bytes"]
            ):
                raise ValueError(
                    f"Retained FASTQ differs from completed mapping record: {path}"
                )
        records.append(
            {
                "sample_id": sample,
                "mapping_metrics_sha256": sha256(metrics_path),
                "mapping_provenance_sha256": sha256(provenance_path),
                "retained_pairs": metrics["retained_pairs"],
                "reads": provenance["outputs"],
            }
        )
    return {
        "status": "PASS",
        "samples": len(rows),
        "minimum_assembly_input_pairs": 0,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["plan", "validate"])
    parser.add_argument(
        "--config", default="config/assembly_completion.json", type=Path
    )
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    if args.stage == "plan":
        completion, cohort = plan(settings)
        write_table(settings["membership"], completion, list(completion[0]))
        write_table(settings["cohort"], cohort, list(cohort[0]))
        print(
            f"{len(completion)} assembly attempts; {len(cohort)} unique cohort libraries"
        )
    else:
        rows = table(settings["membership"])
        report = validate_inputs(rows)
        root = Path(settings["work_root"])
        fields = [
            "sample_id",
            "host_route",
            "strategy",
            "retained_pairs",
            "minimum_assembly_input_pairs",
            "assembly_eligible",
        ]
        write_table(root / "phase3_cohort/input_eligibility.tsv", rows, fields)
        destination = root / "provenance/selected_inputs.json"
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Validated {len(rows)} retained host-depleted inputs")


if __name__ == "__main__":
    main()
