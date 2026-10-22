#!/usr/bin/env python3
"""Compare accepted C. acnes detections with a bounded human-compatible audit.

Counts describe selected FASTQ prefixes. This descriptive comparison neither
tests a contamination mechanism nor converts mapping into human identification.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import statistics
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

VALID = {"validated", "high_confidence"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a small completed reporting table."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def index_rows(rows: list[dict], key: str) -> dict[str, dict]:
    """Reject duplicate or missing keys before joining biological specimens."""
    output = {row[key]: row for row in rows}
    if len(output) != len(rows) or "" in output:
        raise ValueError(f"Nonunique or empty {key}")
    return output


def join_samples(
    metadata: list[dict], bacteria: list[dict], human: list[dict], config: dict
) -> list[dict]:
    """Keep accepted bacterial evidence distinct from mapping audit outcomes."""
    meta = index_rows(metadata, "analysis_sample_id")
    audit = index_rows(human, "sample_id")
    if set(meta) != set(audit) or len(meta) != config["expected_libraries"]:
        raise ValueError("Completed human audit must match the full metadata cohort")
    focal = index_rows(
        [
            row
            for row in bacteria
            if row["target_id"] == config["target_id"]
            and row["grade"] in VALID
            and row["role"] != "decoy"
        ],
        "sample_id",
    )
    if set(focal) - set(meta):
        raise ValueError("Unmatched accepted Cutibacterium specimen")
    output = []
    for sample, row in sorted(meta.items()):
        h = audit[sample]
        if h["selection"] != "equal_per_source_file_prefix":
            raise ValueError("Expected the bounded per-source-prefix human audit")
        raw, qc, matches, unique = (
            int(h[key])
            for key in (
                "raw_selected_pairs",
                "qc_pairs",
                "stringent_human_compatible_pairs",
                "distinct_sequence_pairs",
            )
        )
        if not 0 <= unique <= matches <= qc <= raw or not raw or not qc:
            raise ValueError(f"Invalid selected-prefix counts for {sample}")
        b = focal.get(sample, {})
        output.append(
            {
                "sample_id": sample,
                "library_id": row["library_id"],
                **{
                    key: row[key]
                    for key in (
                        "study",
                        "species_current",
                        "ocean_region",
                        "locality",
                        "tissue",
                    )
                },
                "cutibacterium_accepted": str(bool(b)).lower(),
                "cutibacterium_target": config["target_id"],
                "cutibacterium_grade": b.get("grade", "no_accepted_detection"),
                "cutibacterium_breadth_if_accepted": b.get("breadth", ""),
                "cutibacterium_reads_if_accepted": b.get("reads", ""),
                "cutibacterium_assembly_support_bases_if_accepted": b.get(
                    "assembly_support_bases", ""
                ),
                "raw_prefix_pairs": raw,
                "qc_prefix_pairs": qc,
                "human_compatible_pairs": matches,
                "distinct_human_compatible_sequence_pairs": unique,
                "human_compatible_pair_detected": str(matches > 0).lower(),
                "human_compatible_pairs_per_million_raw_prefix_pairs": matches
                * 1e6
                / raw,
                "human_compatible_pairs_per_million_qc_prefix_pairs": matches
                * 1e6
                / qc,
                "human_reference_sequences": int(h["reference_sequences"]),
                "human_reference_1mb_windows": int(h["occupied_1Mb_windows"]),
                "previously_flagged_specimen": str(
                    sample == config["previously_flagged_sample"]
                ).lower(),
                "human_evidence": "human-compatible mapping; not adjusted by representative nt review",
            }
        )
    return output


def summarize(rows: list[dict], config: dict) -> list[dict]:
    """Report distributions and denominators, without an association test."""
    output = []
    scopes = [
        ("all_specimens", rows),
        (
            "omit_previously_flagged",
            [
                row
                for row in rows
                if row["sample_id"] != config["previously_flagged_sample"]
            ],
        ),
    ]
    studies = ["all_studies"] + sorted({row["study"] for row in rows})
    for scope, scoped in scopes:
        for study in studies:
            selected = [
                row for row in scoped if study == "all_studies" or row["study"] == study
            ]
            for accepted in (True, False):
                subset = [
                    row
                    for row in selected
                    if (row["cutibacterium_accepted"] == "true") == accepted
                ]
                if not subset:
                    continue
                ppm = [
                    row["human_compatible_pairs_per_million_raw_prefix_pairs"]
                    for row in subset
                ]
                pairs = [row["human_compatible_pairs"] for row in subset]
                output.append(
                    {
                        "scope": scope,
                        "study": study,
                        "cutibacterium_evidence": "accepted_detection"
                        if accepted
                        else "no_accepted_detection",
                        "libraries": len(subset),
                        "libraries_with_human_compatible_pairs": sum(
                            value > 0 for value in pairs
                        ),
                        "libraries_without_human_compatible_pairs": sum(
                            value == 0 for value in pairs
                        ),
                        "human_compatible_pairs_median": statistics.median(pairs),
                        "human_compatible_pairs_min": min(pairs),
                        "human_compatible_pairs_max": max(pairs),
                        "human_compatible_pairs_total": sum(pairs),
                        "distinct_sequence_pairs_total": sum(
                            row["distinct_human_compatible_sequence_pairs"]
                            for row in subset
                        ),
                        "raw_prefix_pairs_total": sum(
                            row["raw_prefix_pairs"] for row in subset
                        ),
                        "human_compatible_pairs_per_million_raw_median": statistics.median(
                            ppm
                        ),
                        "human_compatible_pairs_per_million_raw_min": min(ppm),
                        "human_compatible_pairs_per_million_raw_max": max(ppm),
                    }
                )
    return output


def write_table(path: Path, rows: list[dict]) -> None:
    """Write explicit numeric denominators and preserved identifiers."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(config_path: Path) -> None:
    """Join completed tables, export prose values, and preserve input identities."""
    started = datetime.now(timezone.utc)
    config = json.loads(config_path.read_text())
    inputs = {
        key: Path(config[key])
        for key in ("metadata", "bacterial_detections", "human_summary")
    }
    rows = join_samples(*(read_tsv(path) for path in inputs.values()), config)
    summary = summarize(rows, config)
    positives = [row for row in rows if row["cutibacterium_accepted"] == "true"]
    out = Path(config["output_dir"])
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / "sample_comparison.tsv", rows)
    write_table(out / "cutibacterium_specimens.tsv", positives)
    write_table(out / "group_summary.tsv", summary)
    payload = {
        "libraries": len(rows),
        "cutibacterium_positive_libraries": len(positives),
        "cutibacterium_positive_by_study": dict(
            sorted(Counter(row["study"] for row in positives).items())
        ),
        "cutibacterium_positive_by_tissue": dict(
            sorted(Counter(row["tissue"] for row in positives).items())
        ),
        "cutibacterium_positive_by_ocean_region": dict(
            sorted(
                Counter(
                    row["ocean_region"] or "not_recorded_as_ocean_region"
                    for row in positives
                ).items()
            )
        ),
        "group_summary": summary,
        "interpretation": [
            "Descriptive comparison; no formal contamination or association test.",
            "Human-compatible counts describe nonrandom, equally allocated source-file prefixes; they are not whole-library abundance estimates or human-specific confirmations.",
            "No retained human-compatible pair does not establish absence of human DNA.",
            "Cutibacterium status is an accepted bacterial-genome detection; a negative status does not establish absence of the bacterium.",
            "Bacterial contamination can originate from reagents or handling without detectable human DNA; co-occurrence does not establish a shared source.",
            "The known formalin specimen is omitted in a separate sensitivity summary because it was previously flagged by independent background controls.",
        ],
    }
    (out / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    values = {"humanCutibacteriumLibraries": str(len(positives))}
    for scope, prefix in [
        ("all_specimens", "humanCutib"),
        ("omit_previously_flagged", "humanCutibWithoutFlagged"),
    ]:
        for state, suffix in [
            ("accepted_detection", "Detected"),
            ("no_accepted_detection", "NotDetected"),
        ]:
            row = next(
                row
                for row in summary
                if row["scope"] == scope
                and row["study"] == "all_studies"
                and row["cutibacterium_evidence"] == state
            )
            for key, ending in [
                ("libraries", "Libraries"),
                ("libraries_with_human_compatible_pairs", "CompatibleLibraries"),
                ("human_compatible_pairs_median", "MedianCompatiblePairs"),
                ("human_compatible_pairs_max", "MaxCompatiblePairs"),
            ]:
                values[prefix + suffix + ending] = f"{row[key]:g}"
    (out / "values.tex").write_text(
        "% Human-compatible mapping in selected prefixes; descriptive only.\n"
        + "\n".join(
            "\\newcommand{\\" + key + "}{" + value + "}"
            for key, value in values.items()
        )
        + "\n"
    )
    sources = [
        config_path,
        Path("scripts/summarize_human_cutibacterium.py"),
        Path("workflow/human_cutibacterium.smk"),
    ]
    for source in sources:
        destination = out / "provenance/source" / source
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    outputs = [
        out / name
        for name in [
            "sample_comparison.tsv",
            "cutibacterium_specimens.tsv",
            "group_summary.tsv",
            "summary.json",
            "values.tex",
        ]
    ]
    provenance = {
        "analysis_id": "human-compatible-prefix-cutibacterium-comparison",
        "run_id": started.strftime("%Y%m%dT%H%M%SZ"),
        "started_at_utc": started.isoformat(),
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "command": shlex.join([sys.executable, *sys.argv]),
        "python_version": platform.python_version(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "base_git_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "inputs_sha256": {str(path): sha256(path) for path in inputs.values()},
        "source_sha256": {str(path): sha256(path) for path in sources},
        "source_snapshot": str(out / "provenance/source"),
        "outputs_sha256": {str(path): sha256(path) for path in outputs},
        "checks": {
            "unique_full_cohort_join": len(rows),
            "human_count_denominators_reconcile": True,
        },
        "parameters": config,
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(
        json.dumps(
            {
                "libraries": len(rows),
                "cutibacterium_positive_libraries": len(positives),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/human_cutibacterium.json")
    )
    run(parser.parse_args().config)
