#!/usr/bin/env python3
"""Check interpretation-only changes and validate refreshed cohort summaries."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from prepare_coauthor_reporting import digest, read


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    work = Path(config["work_root"])
    output = work / "phase6_analysis/v3"
    baseline = Path(config["baseline_root"])
    record = json.loads((work / "prepared.json").read_text())
    for path, expected in record["inputs_sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Input changed during reporting: {path}")
    for path, expected in record["outputs_sha256"].items():
        if digest(Path(path)) != expected:
            raise ValueError(f"Prepared output changed: {path}")
    originals = read(baseline / "grades/viral_grades.tsv")
    reviewed = read(output / "grades/viral_grades.tsv")
    if len(originals) != len(reviewed):
        raise ValueError("Viral evidence matrix changed size")
    for before, after in zip(originals, reviewed):
        if {k: v for k, v in before.items() if k != "catalog_class"} != {
            k: v for k, v in after.items() if k != "catalog_class"
        }:
            raise ValueError("A viral mapping measurement or grade changed")
    changes = []
    fields = [
        "contamination_tests",
        "nanomia_contamination_tests",
        "bacterial_incidence",
        "viral_incidence",
        "grade_summary",
        "phage_host_links",
    ]
    for threshold in [5, 10, 20]:
        directory = output / (
            "primary" if threshold == 10 else f"sensitivity/breadth_{threshold}pct"
        )
        for name in fields:
            src = output / f"review_tests/{threshold}/{name}.tsv"
            dst = (
                output / "grades/contamination_tests.tsv"
                if threshold == 10 and name == "contamination_tests"
                else directory / f"{name}.tsv"
            )
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            previous = baseline / dst.relative_to(output)
            changes.append(
                {
                    "output": str(dst),
                    "sha256": digest(dst),
                    "baseline_sha256": digest(previous),
                }
            )
    primary = [
        str(output / f"primary/{name}.tsv")
        for name in fields[2:]
        + [
            "physalia_design",
            "physalia_models",
            "physalia_permanova",
            "physalia_region_permutation",
        ]
    ]
    thresholds = [
        str(output / f"sensitivity/breadth_{threshold}pct/{name}.tsv")
        for threshold in [5, 20]
        for name in fields
        + [
            "physalia_design",
            "physalia_models",
            "physalia_permanova",
            "physalia_region_permutation",
        ]
    ]
    sensitivity = [
        str(output / f"sensitivity/{name}.tsv")
        for name in ["threshold_sensitivity", "leave_one_study_out"]
    ]
    command = [
        sys.executable,
        "workflow/scripts/validate_phase6.py",
        "--bacterial",
        str(output / "grades/bacterial_grades.tsv"),
        "--viral",
        str(output / "grades/viral_grades.tsv"),
        "--contamination",
        str(output / "grades/contamination_tests.tsv"),
        "--metadata",
        str(output / "primary/library_flowcells.tsv"),
        "--nanomia",
        str(output / "primary/nanomia_contamination_tests.tsv"),
        "--thresholds",
        *thresholds,
        "--primary",
        *primary,
        "--sensitivity",
        *sensitivity,
        "--samples",
        config["samples"],
        "--output",
        str(work / "validated.done"),
    ]
    subprocess.run(command, check=True)
    (work / "validation.json").write_text(
        json.dumps(
            {
                "status": "PASS",
                "validation_command": command,
                "changed_summaries": changes,
                "mapping_values_and_grades_unchanged": True,
            },
            indent=2,
        )
        + "\n"
    )
    for name in [
        "validate_coauthor_reporting.py",
        "../workflow/scripts/test_contamination_flowcell.py",
        "../workflow/scripts/summarize_incidence.py",
        "../workflow/scripts/phase6_metadata.py",
        "../workflow/scripts/validate_phase6.py",
    ]:
        source = Path("scripts") / name
        shutil.copyfile(source, work / "provenance" / source.name)


if __name__ == "__main__":
    main()
