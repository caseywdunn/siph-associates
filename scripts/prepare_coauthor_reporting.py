#!/usr/bin/env python3
"""Apply interpretation changes to completed mapping evidence, retaining target IDs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def apply_classes(
    rows: list[dict[str, str]], overrides: dict[str, dict[str, str]], key: str
) -> list[dict[str, str]]:
    """Change only interpretation fields; retain all sequence/evidence identities."""
    result = []
    found = set()
    for original in rows:
        row = dict(original)
        target = row[key]
        if target in overrides:
            change = overrides[target]
            if row["catalog_class"] != change["original_class"]:
                raise ValueError(f"Unexpected existing class for {target}")
            if "phylum" in row and row["phylum"] != change["expected_phylum"]:
                raise ValueError(f"Unexpected taxonomy for {target}")
            row["catalog_class"] = change["reporting_class"]
            if "endogenous_reasons" in row:
                row["endogenous_reasons"] = ";".join(
                    filter(
                        None, [row["endogenous_reasons"], "RNA_virus_in_DNA_library"]
                    )
                )
            found.add(target)
        result.append(row)
    if found != set(overrides):
        raise ValueError(f"Missing override targets: {set(overrides) - found}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    baseline = Path(config["baseline_root"])
    work = Path(config["work_root"])
    output = work / "phase6_analysis/v3"
    source = Path(config["catalog_source"])
    catalog = Path(config["catalog_output"])
    overrides_rows = read(Path(config["classifications"]))
    overrides = {r["target_id"]: r for r in overrides_rows}
    if len(overrides) != len(overrides_rows):
        raise ValueError("Duplicate viral classification overrides")
    inputs = {
        str(p): digest(p)
        for p in [
            args.config,
            Path(config["classifications"]),
            Path(config["baseline_validation"]),
        ]
    }
    outputs = {}
    reused = []

    def checked_copy(src: Path, dst: Path) -> None:
        inputs[str(src)] = digest(src)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        outputs[str(dst)] = digest(dst)
        if outputs[str(dst)] != inputs[str(src)]:
            raise ValueError(f"Copy differs: {dst}")
        reused.append(
            {"source": str(src), "output": str(dst), "sha256": inputs[str(src)]}
        )

    old_flow = baseline / "primary/library_flowcells.tsv"
    new_flow = output / "primary/library_flowcells.tsv"
    for path in [old_flow, new_flow]:
        inputs[str(path)] = digest(path)

    def physalia(records: list[dict[str, str]]) -> list[dict[str, str]]:
        return [r for r in records if r["host_species"].startswith("Physalia ")]

    if physalia(read(old_flow)) != physalia(read(new_flow)):
        raise ValueError("Physalia metadata changed; regional model reuse is invalid")
    checked_copy(
        baseline / "grades/bacterial_grades.tsv", output / "grades/bacterial_grades.tsv"
    )
    for suffix in ["primary", "sensitivity/breadth_5pct", "sensitivity/breadth_20pct"]:
        for name in [
            "physalia_design",
            "physalia_models",
            "physalia_permanova",
            "physalia_region_permutation",
        ]:
            checked_copy(
                baseline / suffix / f"{name}.tsv", output / suffix / f"{name}.tsv"
            )
    for name in ["threshold_sensitivity", "leave_one_study_out"]:
        checked_copy(
            baseline / "sensitivity" / f"{name}.tsv",
            output / "sensitivity" / f"{name}.tsv",
        )
    for name in [
        "bacterial_catalog.manifest.tsv",
        "mag_catalog.tsv",
        "host_links.tsv",
        "catalog_id_crosswalk.tsv",
    ]:
        checked_copy(source / name, catalog / name)
    manifests = [
        source / f"viral_{kind}.manifest.tsv"
        for kind in ["associate", "endogenous_candidate"]
    ]
    original_manifests = [r for p in manifests for r in read(p)]
    for path in manifests:
        inputs[str(path)] = digest(path)
    updated_manifests = apply_classes(original_manifests, overrides, "votu_id")
    if len({r["votu_id"] for r in updated_manifests}) != len(updated_manifests):
        raise ValueError("Viral target identities are not unique")
    for kind in ["associate", "endogenous_candidate"]:
        write(
            catalog / f"viral_{kind}.manifest.tsv",
            [r for r in updated_manifests if r["catalog_class"] == kind],
        )
    grades_path = baseline / "grades/viral_grades.tsv"
    inputs[str(grades_path)] = digest(grades_path)
    original_grades = read(grades_path)
    updated_grades = apply_classes(original_grades, overrides, "target_id")
    write(output / "grades/viral_grades.tsv", updated_grades)
    removed = [
        r
        for r in original_grades
        if r["target_id"] in overrides
        and r["grade"] in {"validated", "high_confidence"}
    ]
    write(catalog / "reclassified_detections.tsv", removed)
    summary = {
        "catalog_counts": dict(Counter(r["catalog_class"] for r in updated_manifests)),
        "reclassified_targets": sorted(overrides),
        "reclassified_accepted_detections": len(removed),
        "reclassified_assembly_supported": sum(
            r["grade"] == "high_confidence" for r in removed
        ),
        "target_ids_and_sequence_checksums_unchanged": True,
        "mapping_inputs_and_grades_unchanged": True,
        "physalia_models_reused": "Physalia metadata and bacterial evidence unchanged; original validated run retained.",
        "integration_status": "Possible endogenous candidates; genomic integration is not demonstrated.",
    }
    (catalog / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    for path in [
        output / "grades/viral_grades.tsv",
        *catalog.glob("*.tsv"),
        catalog / "summary.json",
    ]:
        outputs[str(path)] = digest(path)
    code = Path(__file__)
    record = {
        "summary": summary,
        "command": sys.argv,
        "python": sys.version,
        "source_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "code_sha256": digest(code),
        "inputs_sha256": inputs,
        "outputs_sha256": outputs,
        "reused_outputs": reused,
    }
    (work / "provenance").mkdir(parents=True, exist_ok=True)
    shutil.copyfile(code, work / "provenance/prepare_coauthor_reporting.py")
    for path in [work / "prepared.json", catalog / "provenance.json"]:
        path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
