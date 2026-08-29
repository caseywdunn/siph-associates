#!/usr/bin/env python3
"""Validate Phase-2 design and optionally the pilot or complete cohort outputs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "data" / "results"
REPORT = ROOT / "data" / "metadata" / "phase2_validation.txt"

parser = argparse.ArgumentParser()
parser.add_argument("--scope", choices=("static", "pilot", "cohort"), default="static")
args = parser.parse_args()
errors = []


def require(condition, message):
    if not condition:
        errors.append(message)


def table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


samples = table(ROOT / "config" / "samples.tsv")
sample_by_id = {row["sample_id"]: row for row in samples if row["include_primary"].lower() == "true"}
pilot = table(ROOT / "config" / "phase2_pilot.tsv")
exclusions = table(ROOT / "config" / "phase2_exclusions.tsv")
pilot_ids = [row["sample_id"] for row in pilot]
excluded_ids = {row["sample_id"] for row in exclusions}

require(len(sample_by_id) == 205, f"expected 205 included samples, found {len(sample_by_id)}")
require(6 <= len(pilot_ids) <= 10, f"pilot must contain 6--10 samples, found {len(pilot_ids)}")
require(len(pilot_ids) == len(set(pilot_ids)), "pilot sample IDs are duplicated")
require(set(pilot_ids).issubset(sample_by_id), "pilot contains IDs absent from production sample table")
require({sample_by_id[s]["study"] for s in pilot_ids} == {"Church2025", "Ahuja2024", "Ahuja2026"},
        "pilot does not cover all three studies")
require({sample_by_id[s]["host_route"] for s in pilot_ids} == {"P_physalis", "N_septata", "none"},
        "pilot does not cover all three host routes")
pilot_depths = [int(sample_by_id[s]["read_pairs"]) for s in pilot_ids]
require(any(depth > 400_000_000 for depth in pilot_depths), "pilot lacks a capped library")
require(any(depth < 400_000_000 for depth in pilot_depths), "pilot lacks an uncapped library")
require(excluded_ids.issubset(sample_by_id), "exclusion table contains an unknown sample ID")
require(all(row.get("reason", "").strip() for row in exclusions), "an exclusion lacks a reason")

snake_text = "\n".join(path.read_text() for path in [ROOT / "Snakefile", *sorted((ROOT / "workflow" / "rules").glob("*.smk"))])
for rule in ("validate_screen_sample", "aggregate_phase2", "screen_pilot", "screen_cohort"):
    require(f"rule {rule}:" in snake_text, f"missing Phase-2 rule: {rule}")

scope_ids = []
if args.scope == "pilot":
    scope_ids = pilot_ids
elif args.scope == "cohort":
    scope_ids = [sample for sample in sorted(sample_by_id) if sample not in excluded_ids]

validation_records = []
if scope_ids:
    for sample in scope_ids:
        path = WORK / "validation" / "screens" / f"{sample}.json"
        if not path.is_file():
            errors.append(f"missing screen validation: {sample}")
            continue
        try:
            record = json.loads(path.read_text())
            validation_records.append(record)
            require(record.get("status") == "PASS", f"screen validation is not PASS: {sample}")
            require(record.get("sample_id") == sample, f"screen validation sample mismatch: {sample}")
            for field in ("library_id", "specimen_id", "study", "raw_pairs", "pairs_after_fastp", "screen_rows"):
                require(field in record, f"screen validation lacks {field}: {sample}")
        except Exception as exc:
            errors.append(f"cannot parse screen validation {sample}: {exc}")
        required_paths = [
            WORK / "trim" / f"{sample}.fastp.json",
            WORK / "provenance" / "trim" / f"{sample}.json",
            WORK / "screens" / "kraken" / f"{sample}.report.tsv",
            WORK / "screens" / "kraken" / f"{sample}.assignments.tsv.gz",
            WORK / "screens" / "kraken" / f"{sample}.classified_1.fastq.gz",
            WORK / "screens" / "kraken" / f"{sample}.classified_2.fastq.gz",
            WORK / "screens" / "bracken" / f"{sample}.G.tsv",
            WORK / "screens" / "bracken" / f"{sample}.S.tsv",
            WORK / "screens" / "sylph" / f"{sample}.profile.tsv",
            WORK / "screens" / "phyloflash" / f"{sample}.tar.gz",
        ]
        for path in required_paths:
            require(path.is_file() and path.stat().st_size > 0, f"missing or empty output: {path.relative_to(ROOT)}")

    aggregate_scope = "pilot" if args.scope == "pilot" else "cohort"
    qc_path = WORK / "aggregation" / aggregate_scope / "library_qc.tsv"
    nomination_path = WORK / "aggregation" / aggregate_scope / "candidate_nominations.tsv"
    stage_path = WORK / "stages" / f"screen_{aggregate_scope}.done"
    require(qc_path.is_file(), f"missing {aggregate_scope} QC aggregation")
    require(nomination_path.is_file(), f"missing {aggregate_scope} nomination aggregation")
    require(stage_path.is_file(), f"missing screen_{aggregate_scope} sentinel")
    if qc_path.is_file():
        qc = table(qc_path)
        require(len(qc) == len(scope_ids), f"{aggregate_scope} QC row count mismatch")
        require({row.get("sample_id") for row in qc} == set(scope_ids), f"{aggregate_scope} QC sample set mismatch")
    if nomination_path.is_file():
        nominations = table(nomination_path)
        require(bool(nominations), f"{aggregate_scope} nomination table has no data rows")
        if nominations:
            require({"sample_id", "study", "source", "rank", "candidate_id", "candidate_name"}.issubset(nominations[0]),
                    f"{aggregate_scope} nomination table lacks required columns")
    if stage_path.is_file():
        stage = dict(line.split("\t", 1) for line in stage_path.read_text().splitlines() if "\t" in line)
        require(stage.get("status") == "PASS", f"screen_{aggregate_scope} sentinel is not PASS")
        require(stage.get("samples") == str(len(scope_ids)), f"screen_{aggregate_scope} sample count mismatch")

metrics = {
    "status": "PASS" if not errors else "FAIL",
    "scope": args.scope,
    "production_samples": len(sample_by_id),
    "pilot_samples": len(pilot_ids),
    "excluded_samples": len(excluded_ids),
    "validated_samples": len(validation_records),
    "errors": len(errors),
}
REPORT.parent.mkdir(parents=True, exist_ok=True)
REPORT.write_text("".join(f"{key}\t{value}\n" for key, value in metrics.items()) +
                  "".join(f"ERROR\t{error}\n" for error in errors))
print(REPORT.read_text(), end="")
if errors:
    raise SystemExit(1)
