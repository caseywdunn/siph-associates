#!/usr/bin/env python3
"""Report host-depletion counts from completed mapping summaries, without FASTQ I/O."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

FIELDS = ["species", "sample", "reads_ingested", "percent_removed", "reads_remaining"]
SOURCE_FILES = [
    "scripts/summarize_host_depletion.py",
    "workflow/rules/host_depletion_report.smk",
    "workflow/host_depletion_report.smk",
]


def sha256(path: str | Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_table(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def selected_samples(settings: dict) -> dict[str, dict]:
    """Select reference-bearing libraries; exclude other assembly-input routes."""
    selected = {}
    for source in settings["sources"]:
        for row in read_table(source["samples"]):
            if row["include_primary"].lower() != "true":
                continue
            route = row["host_route"]
            if route in ("", "none"):
                continue
            sample = row["sample_id"]
            if sample in selected:
                raise ValueError(f"Duplicate reference-bearing sample: {sample}")
            root = Path(source["results"])
            selected[sample] = {
                "host_route": route,
                "metrics": str(root / "phase3_cohort/inputs" / f"{sample}.json"),
                "mapping": str(
                    root / "provenance/phase3_cohort/inputs" / f"{sample}.json"
                ),
                "trim": str(root / "provenance/trim" / f"{sample}.json"),
            }
    if not selected:
        raise ValueError("No reference-bearing libraries selected")
    unknown = set(settings.get("allowed_pending_samples", [])) - selected.keys()
    if unknown:
        raise ValueError(f"Unknown pending samples: {sorted(unknown)}")
    return selected


def report_inputs(settings: dict, require_complete: bool = False) -> list[str]:
    """Declare all small input files, allowing only explicitly pending mappings."""
    paths = [settings["metadata"], *(s["samples"] for s in settings["sources"])]
    pending = set() if require_complete else set(settings["allowed_pending_samples"])
    for sample, source in selected_samples(settings).items():
        paths.append(source["trim"])
        for key in ("metrics", "mapping"):
            path = source[key]
            if sample not in pending or Path(path).is_file():
                paths.append(path)
    return paths


def integer_count(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} must be a nonnegative integer")
    return value


def build_rows(
    settings: dict, require_complete: bool = False
) -> tuple[list[dict], list[str], dict[str, str]]:
    """Join current names to stable analysis IDs and validate mapping provenance."""
    checksums = {}

    def read_json(path: str) -> dict:
        data = Path(path).read_bytes()
        checksums[path] = hashlib.sha256(data).hexdigest()
        return json.loads(data)

    metadata_rows = read_table(settings["metadata"])
    metadata = {r["analysis_sample_id"]: r for r in metadata_rows}
    if len(metadata) != len(metadata_rows):
        raise ValueError("Duplicate analysis_sample_id in metadata")
    for path in [settings["metadata"], *(s["samples"] for s in settings["sources"])]:
        checksums[path] = sha256(path)
    allowed_pending = (
        set() if require_complete else set(settings["allowed_pending_samples"])
    )
    rows, pending = [], []
    for sample, source in selected_samples(settings).items():
        if sample not in metadata or not metadata[sample]["species_current"].strip():
            raise ValueError(f"Missing species metadata for {sample}")
        trim = read_json(source["trim"])
        if trim["sample_id"] != sample:
            raise ValueError(f"Trim sample mismatch for {sample}")
        input_pairs = integer_count(trim["pairs_after_fastp"], f"{sample} input_pairs")
        row = {
            "species": metadata[sample]["species_current"],
            "sample": sample,
            "reads_ingested": input_pairs * 2,
            "percent_removed": "NA",
            "reads_remaining": "NA",
        }
        # Both records are written after successful mapping and retained-read counting.
        # A pilot or a lone intermediate record must not supply final counts.
        completed = all(Path(source[k]).is_file() for k in ("metrics", "mapping"))
        if not completed:
            if sample not in allowed_pending:
                raise FileNotFoundError(f"Missing completed host mapping for {sample}")
            for key in ("metrics", "mapping"):
                if Path(source[key]).is_file():
                    read_json(source[key])
            pending.append(sample)
            rows.append(row)
            continue
        metrics = read_json(source["metrics"])
        mapping = read_json(source["mapping"])
        for record in (metrics, mapping):
            expected = {
                "sample_id": sample,
                "host_route": source["host_route"],
                "strategy": "reference_depleted",
            }
            if any(record.get(k) != v for k, v in expected.items()):
                raise ValueError(f"Mapping identity/strategy mismatch for {sample}")
        if mapping["upstream_provenance_sha256"] != checksums[source["trim"]]:
            raise ValueError(f"Mapping/trim provenance mismatch for {sample}")
        measured_input = integer_count(metrics["input_pairs"], f"{sample} input_pairs")
        if measured_input != input_pairs:
            raise ValueError(f"Mapping/trim counts disagree for {sample}")
        retained = integer_count(metrics["retained_pairs"], f"{sample} retained_pairs")
        if input_pairs == 0 or retained > input_pairs:
            raise ValueError(f"Invalid mapping count bounds for {sample}")
        if not math.isclose(
            metrics["retained_fraction"],
            retained / input_pairs,
            rel_tol=1e-12,
            abs_tol=1e-12,
        ):
            raise ValueError(f"Retained fraction disagrees with counts for {sample}")
        row["percent_removed"] = f"{100 * (input_pairs - retained) / input_pairs:.4f}"
        row["reads_remaining"] = retained * 2
        rows.append(row)
    rows.sort(key=lambda r: (r["species"], r["sample"]))
    return rows, sorted(pending), checksums


def write_markdown(path: Path, rows: list[dict], pending: list[str]) -> None:
    lines = [
        "# Host-depletion mapping statistics",
        "",
        f"{len(rows)} samples; {len(rows) - len(pending)} completed mappings.",
        "Counts are individual reads (two per pair) after quality filtering.",
        "A pair is retained only when both mates are unmapped; percent removed",
        "counts both reads of every removed pair.",
        "",
    ]
    if pending:
        lines.extend(
            [
                "**Pending final mapping summaries:** " + ", ".join(pending) + ".",
                "`NA` denotes unavailable counts, not zero. Input counts are from",
                "completed trimming; pilot counts are not substituted.",
                "",
            ]
        )
    lines.extend(
        [
            "| Sample species | Sample | Reads ingested | Percent removed | Reads remaining |",
            "| --- | --- | ---: | ---: | ---: |",
        ]
    )
    for row in rows:
        cells = [
            f"{row[k]:,}" if isinstance(row[k], int) else str(row[k]) for k in FIELDS
        ]
        lines.append("| " + " | ".join(cells) + " |")
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    rows, pending, checksums = build_rows(settings, args.require_complete)
    destination = Path(settings["output_dir"])
    destination.mkdir(parents=True, exist_ok=True)
    table = destination / "mapping_stats.tsv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    markdown = destination / "mapping_stats.md"
    write_markdown(markdown, rows, pending)
    # Preserve the small, possibly uncommitted producer alongside its run record.
    for source in [*SOURCE_FILES, str(args.config)]:
        snapshot = destination / "provenance/source" / source
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, snapshot)
        checksums[source] = sha256(snapshot)
    provenance = {
        "generated_utc": datetime.now(UTC).isoformat(),
        "command": [sys.executable, *sys.argv],
        "python": platform.python_version(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "status": "partial" if pending else "complete",
        "samples": len(rows),
        "completed_mappings": len(rows) - len(pending),
        "pending_samples": pending,
        "require_complete": args.require_complete,
        "units": "individual reads; two reads per pair",
        "reads_ingested": "quality-filtered reads entering host depletion",
        "percent_removed": "100 * (input_pairs - retained_pairs) / input_pairs",
        "retention": "both mates unmapped (samtools fastq -f 12)",
        "sample_identifier": "analysis_sample_id from curated metadata",
        "inputs_sha256": checksums,
        "outputs_sha256": {str(p): sha256(p) for p in (table, markdown)},
    }
    (destination / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(
        f"{table}: {len(rows)} samples, {len(rows) - len(pending)} complete, "
        f"{len(pending)} pending"
    )


if __name__ == "__main__":
    main()
