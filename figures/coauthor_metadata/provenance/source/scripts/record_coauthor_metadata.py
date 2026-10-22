#!/usr/bin/env python3
"""Preserve the completed metadata analysis and its small source/input identities."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from analyze_collection_context import INPUTS, read_table


def digest(path: Path) -> str:
    """Hash the actual file supporting this execution."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_values(out: Path, tests: list[dict]) -> None:
    """Expose reported numbers without asking the manuscript to recompute them."""
    values = {}
    for contrast, prefix in [
        ("utriculus_pacific_regions", "metadataWithinHost"),
        ("sw_pacific_host_species", "metadataWithinRegion"),
    ]:
        for grade, suffix in [
            ("grade", ""),
            ("grade_at_5pct", "FivePct"),
            ("grade_at_20pct", "TwentyPct"),
        ]:
            row = next(
                row
                for row in tests
                if row["contrast"] == contrast and row["grade_column"] == grade
            )
            if row["status"] != "tested":
                raise ValueError(
                    f"Cannot export a manuscript statistic for {contrast}/{grade}"
                )
            values[prefix + suffix + "Libraries"] = row["included_libraries"]
            values[prefix + suffix + "RtwoPercent"] = (
                f"{100 * float(row['r_squared']):.1f}"
            )
            values[prefix + suffix + "P"] = f"{float(row['p_value']):.4f}"
            values[prefix + suffix + "Q"] = f"{float(row['q_value']):.4f}"
    tissue = read_table(out / "tissue_incidence.tsv")
    for category, prefix in [
        ("zooid", "metadataFishZooid"),
        ("tentacle", "metadataFishTentacle"),
    ]:
        row = next(
            row
            for row in tissue
            if row["host_group"] == "Physalia" and row["tissue_group"] == category
        )
        values[prefix + "Positive"] = row["fish_12s_positive_libraries"]
        values[prefix + "Libraries"] = row["total_libraries"]
        values[prefix + "Pairs"] = row["fish_12s_matching_pairs"]
    descriptive = json.loads((out / "descriptive_values.json").read_text())
    values["metadataFishAllPairs"] = str(descriptive["fish_12s_matching_pairs"])
    values["metadataDepthLibraries"] = str(descriptive["numeric_depths"])
    for host, prefix in [
        ("Physalia", "metadataReferencePhysalia"),
        ("Nanomia septata", "metadataReferenceNanomia"),
        ("Halistemma rubrum", "metadataReferenceHalistemma"),
    ]:
        row = next(
            row
            for row in read_table(out / "host_depletion_summary.tsv")
            if row["scope"] == "reference_group" and row["host"] == host
        )
        values[prefix + "Libraries"] = row["libraries"]
        for key, suffix in [
            ("percent_removed_min", "MinPercentRemoved"),
            ("percent_removed_median", "MedianPercentRemoved"),
            ("percent_removed_max", "MaxPercentRemoved"),
        ]:
            values[prefix + suffix] = f"{float(row[key]):.2f}"
    for row in read_table(out / "historical_host_concordance.tsv"):
        if row["host_group"] != "all":
            continue
        prefix = "metadataHistoricalHost" + row["domain"].capitalize()
        for key, suffix in [
            ("libraries", "Libraries"),
            ("present_host_depleted", "Depleted"),
            ("present_full_reads", "Full"),
            ("both", "Shared"),
            ("only_host_depleted", "OnlyDepleted"),
            ("only_full_reads", "OnlyFull"),
        ]:
            values[prefix + suffix] = row[key]
    text = "% Generated from the completed coauthor metadata analysis.\n"
    text += "% HistoricalHost macros describe fixed v2 controls; viral controls include endogenous candidates.\n"
    text += "\n".join(
        "\\newcommand{\\" + name + "}{" + value + "}" for name, value in values.items()
    )
    (out / "values.tex").write_text(text + "\n")


def record(config_path: Path) -> None:
    """Validate complete identifiers and freeze source with output hashes."""
    config = json.loads(config_path.read_text())
    out = Path(config["output_dir"])
    context = read_table(out / "specimen_context.tsv")
    ids = [row["sample_id"] for row in context]
    if len(ids) != config["expected_libraries"] or len(ids) != len(set(ids)):
        raise ValueError("Incomplete or duplicated specimen cohort")
    tests = read_table(out / "community_tests.tsv")
    if len(tests) != len(config["grade_columns"]) * len(config["contrasts"]):
        raise ValueError("Missing contrast or threshold result")
    for row in tests:
        if int(row["included_libraries"]) + int(row["excluded_empty_profiles"]) != int(
            row["eligible_libraries"]
        ):
            raise ValueError("Community inclusion counts do not reconcile")
        if row["status"] == "tested" and not 0 <= float(row["p_value"]) <= 1:
            raise ValueError("Invalid tested p-value")
        if row["status"] == "tested" and float(row["residual_df"]) < 2:
            raise ValueError("Tested contrast lacks residual degrees of freedom")
    current_metadata = {
        row["analysis_sample_id"]: row for row in read_table(Path(INPUTS["manifest"]))
    }
    for row in context:
        current = current_metadata[row["sample_id"]]
        if any(
            row[field] != current[field]
            for field in ("species_current", "depth_m", "library_id")
        ):
            raise ValueError(
                "Publication metadata changed during execution; rerun the bounded analysis"
            )
    write_values(out, tests)
    sources = [
        config_path,
        Path("scripts/analyze_coauthor_metadata.py"),
        Path("scripts/analyze_collection_context.py"),
        Path("workflow/scripts/analyze_coauthor_communities.R"),
        Path(__file__).relative_to(Path.cwd()),
        Path("workflow/coauthor_metadata.smk"),
        Path("workflow/rules/coauthor_metadata.smk"),
        Path("envs/stats.yaml"),
    ]
    inputs = list(
        dict.fromkeys(
            [
                *INPUTS.values(),
                config["mapping_summary"],
                config["bacterial_grades"],
                config["fish_summary"],
                config["host_depletion_stats"],
                config["host_depletion_provenance"],
                config["historical_host_concordance"],
            ]
        )
    )
    outputs = sorted(
        path
        for path in out.iterdir()
        if path.is_file() and path.name != "provenance.json"
    )
    if len(outputs) != 23:
        raise ValueError(f"Expected 23 analysis artifacts; found {len(outputs)}")
    for path in sources:
        target = out / "provenance/source" / path
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    provenance = {
        "analysis_id": "coauthor-metadata-patterns",
        "run_id": os.environ.get("SLURM_JOB_ID", "local")
        + "-"
        + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "status": "complete",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": shlex.join([sys.executable, *sys.argv]),
        "workflow": "snakemake -s workflow/coauthor_metadata.smk all --cores 1 --rerun-triggers mtime",
        "source_base_git_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "source_snapshot": str(out / "provenance/source"),
        "source_sha256": {str(path): digest(path) for path in sources},
        "inputs_sha256": {path: digest(Path(path)) for path in inputs},
        "outputs_sha256": {str(path): digest(path) for path in outputs},
        "environment": {
            "python": platform.python_version(),
            "r_session": str(out / "session_info.txt"),
        },
        "parameters": config,
        "checks": {
            "unique_complete_cohort": len(ids),
            "all_contrast_thresholds": len(tests),
            "inclusion_counts_reconcile": True,
        },
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/coauthor_metadata_analysis.json")
    )
    record(parser.parse_args().config)
