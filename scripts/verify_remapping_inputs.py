#!/usr/bin/env python3
"""Verify explicit mapping authorization and reuse the reviewed read inventory."""

from __future__ import annotations

import argparse
import json
import shutil
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from assembly_completion import sha256, table, unique, write_table


def verify_authorization(settings: dict) -> dict:
    authorization = json.loads(Path(settings["authorization"]).read_text())
    if authorization["status"] != "authorized" or authorization["automatic_remapping"]:
        raise ValueError("This remapping round requires explicit user authorization")
    approved = authorization["authorization"]
    if approved["catalog_version"] != settings["catalog_version"]:
        raise ValueError("Authorization does not cover this catalog version")
    for key in (
        "catalog_readiness",
        "catalog_inventory",
        "mapping_read_inventory",
        "mapping_policy",
    ):
        if sha256(settings[key]) != approved[f"{key}_sha256"]:
            raise ValueError(f"Authorized input has changed: {key}")
    ready = json.loads(Path(settings["catalog_readiness"]).read_text())
    if ready["status"] != "PASS" or int(ready["errors"]) != 0:
        raise ValueError("Catalog validation has not passed")
    if (
        Path(ready["catalog_root"]).resolve()
        != Path(settings["mapping_catalog_root"]).resolve()
    ):
        raise ValueError(
            "Mapping points to a different catalog than the readiness report"
        )
    return approved


def verify_read_record(row: dict, record: dict) -> None:
    """Check producer identity, pair accounting and sizes without decompressing reads."""
    sample = row["sample_id"]
    if record["sample_id"] != sample or len(record["outputs"]) != 2:
        raise ValueError(f"Read producer identity mismatch: {sample}")
    if row["host_route"] == "none":
        if int(record["pairs_after_fastp"]) != int(row["read_pairs"]):
            raise ValueError(f"Full trimmed pair count mismatch: {sample}")
        if record["parameters"]["cap_pairs"] != int(row["cap_pairs"]) or record[
            "parameters"
        ]["cap_seed"] != int(row["cap_seed"]):
            raise ValueError(f"Read capping differs from the frozen policy: {sample}")
    elif (
        record.get("strategy") != "reference_depleted"
        or record.get("host_route") != row["host_route"]
    ):
        raise ValueError(f"Host-depleted read source mismatch: {sample}")
    for mate, output in enumerate(record["outputs"], 1):
        path = Path(row[f"r{mate}"])
        if (
            path.name != Path(output["path"]).name
            or path.stat().st_size != output["size_bytes"]
        ):
            raise ValueError(f"Read file differs from its producer record: {path}")


def verify_regeneration(row: dict, baseline: dict, current: dict) -> None:
    """Reject changed raw inputs, capped subsets, or trimming before mapping them."""
    verify_read_record(row, current)
    if baseline["sample_id"] != current["sample_id"]:
        raise ValueError("Regenerated sample identity mismatch")
    for key in (
        "software",
        "raw_pairs_expected",
        "pairs_entering_fastp",
        "pairs_after_fastp",
        "parameters",
    ):
        if baseline[key] != current[key]:
            raise ValueError(
                f"Regenerated trimming differs from its baseline: {row['sample_id']} {key}"
            )
    before = [(r["path"], r["size_bytes"], r["sha256"]) for r in baseline["inputs"]]
    after = [(r["path"], r["size_bytes"], r["sha256"]) for r in current["inputs"]]
    if before != after:
        raise ValueError(f"Raw input identities changed: {row['sample_id']}")
    # fastp/BBTools gzip headers need not be byte-identical across executions.
    # The fixed software, raw hashes, ordered lanes, cap/seed and paired count
    # establish reuse of the same trimming method; new output hashes are retained.


def preflight(settings: dict) -> None:
    approved = verify_authorization(settings)
    rows = table(settings["mapping_read_inventory"])
    inputs = unique(rows, "sample_id")
    cohort = unique(table(settings["samples"]), "sample_id")
    if set(inputs) != set(cohort) or len(inputs) != settings["expected_samples"]:
        raise ValueError("Read inventory must match the complete expanded cohort")
    for record in table(settings["catalog_inventory"]):
        path = Path(record["path"])
        if (
            path.stat().st_size != int(record["bytes"])
            or sha256(path) != record["sha256"]
        ):
            raise ValueError(f"Catalog artifact changed after validation: {path}")
    counts = Counter()
    eligibility, qc, sources = [], [], []
    for row in rows:
        sample = row["sample_id"]
        baseline = Path(row["baseline_provenance"])
        if sha256(baseline) != row["baseline_provenance_sha256"]:
            raise ValueError(f"Read provenance changed: {sample}")
        record = json.loads(baseline.read_text())
        if row["host_route"] != cohort[sample]["host_route"]:
            raise ValueError(f"Cohort host route mismatch: {sample}")
        if (
            int(row["cap_pairs"]) != settings["parameters"]["cap_pairs"]
            or int(row["cap_seed"]) != settings["parameters"]["cap_seed"]
        ):
            raise ValueError(f"Unexpected cap/seed: {sample}")
        if row["availability"] == "ready":
            verify_read_record(row, record)
            if row["host_route"] == "none":
                counts["retained_trimmed"] += 1
                destination = (
                    Path(settings["work_root"])
                    / "provenance/phase5_mapping/trim"
                    / f"{sample}.json"
                )
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(baseline, destination)
            else:
                counts["host_depleted"] += 1
                metrics = (
                    Path(cohort[sample]["input_root"])
                    / "phase3_cohort/inputs"
                    / f"{sample}.json"
                )
                measured = json.loads(metrics.read_text())
                if (
                    measured["sample_id"] != sample
                    or measured["host_route"] != row["host_route"]
                    or int(measured["retained_pairs"]) != int(row["read_pairs"])
                ):
                    raise ValueError(
                        f"Host-depletion pair accounting mismatch: {sample}"
                    )
        else:
            if row["host_route"] != "none" or row["r1"] or row["r2"]:
                raise ValueError(
                    f"Only missing full trimmed reads may be regenerated: {sample}"
                )
            counts["regenerated"] += 1
            paths = row["raw_r1_paths"].split(";") + row["raw_r2_paths"].split(";")
            if paths != [r["path"] for r in record["inputs"]]:
                raise ValueError(f"Raw lane order differs from baseline: {sample}")
            if any(
                Path(r["path"]).stat().st_size != r["size_bytes"]
                for r in record["inputs"]
            ):
                raise ValueError(f"Raw input size changed: {sample}")
            if record["sample_id"] != sample or record["pairs_after_fastp"] != int(
                row["read_pairs"]
            ):
                raise ValueError(f"Regeneration baseline mismatch: {sample}")
        eligibility.append({"sample_id": sample, "retained_pairs": row["read_pairs"]})
        qc.append(
            {
                "sample_id": sample,
                "pairs_after_fastp": record.get("pairs_after_fastp", ""),
            }
        )
        sources.append(
            {
                "sample_id": sample,
                "provenance_sha256": row["baseline_provenance_sha256"],
                "read_pairs": int(row["read_pairs"]),
                "availability": row["availability"],
            }
        )
    for kind in ("regenerated", "host_depleted", "retained_trimmed"):
        if counts[kind] != settings[f"expected_{kind}_libraries"]:
            raise ValueError(f"Unexpected number of {kind} libraries: {counts[kind]}")
    write_table(settings["mapping_eligibility"], eligibility, list(eligibility[0]))
    write_table(settings["mapping_library_qc"], qc, list(qc[0]))
    report = {
        "status": "PASS",
        "checked_utc": datetime.now(UTC).isoformat(),
        "authorization": approved,
        "libraries": len(rows),
        "read_sources": dict(counts),
        "catalog_checksums_verified": True,
        "read_verification": "completed producer records and current sizes; regenerated raw hashes rechecked before mapping",
        "sources": sources,
    }
    Path(settings["mapping_preflight"]).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "sources"}, indent=2))


def regenerated(settings: dict, sample: str) -> None:
    row = unique(table(settings["mapping_read_inventory"]), "sample_id")[sample]
    root = Path(settings["work_root"])
    for mate in (1, 2):
        row[f"r{mate}"] = str(
            root / "phase5_mapping/trimmed" / f"{sample}_R{mate}.fastq.gz"
        )
    provenance = root / "provenance/phase5_mapping/trim" / f"{sample}.json"
    current = json.loads(provenance.read_text())
    baseline = json.loads(Path(row["baseline_provenance"]).read_text())
    verify_regeneration(row, baseline, current)
    report = {
        "status": "PASS",
        "sample_id": sample,
        "pairs": current["pairs_after_fastp"],
        "provenance_sha256": sha256(provenance),
        "outputs": current["outputs"],
    }
    output = Path(settings["mapping_read_checks"]) / f"{sample}.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["preflight", "regenerated"])
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--sample")
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    if args.stage == "preflight":
        preflight(settings)
    else:
        if not args.sample:
            parser.error("regenerated requires --sample")
        regenerated(settings, args.sample)


if __name__ == "__main__":
    main()
