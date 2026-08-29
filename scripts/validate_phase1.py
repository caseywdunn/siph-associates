#!/usr/bin/env python3
"""Validate Phase-1 structure and, optionally, completed smoke outputs."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "data" / "metadata" / "phase1_validation.txt"

parser = argparse.ArgumentParser()
parser.add_argument("--require-smoke", action="store_true")
args = parser.parse_args()
errors = []


def require(condition, message):
    if not condition:
        errors.append(message)


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


production = rows(ROOT / "config" / "samples.tsv")
fixture = rows(ROOT / "tests" / "fixtures" / "samples.tsv")
require(len(production) == 205, f"production samples: expected 205, found {len(production)}")
require(len({row["sample_id"] for row in production}) == len(production), "production sample IDs are duplicated")
require(len(fixture) == 3, f"fixture samples: expected 3, found {len(fixture)}")
require({row["study"] for row in fixture} == {"Church2025", "Ahuja2024", "Ahuja2026"},
        "fixture panel does not cover three studies")
require({row["host_route"] for row in fixture} == {"P_physalis", "N_septata", "none"},
        "fixture panel does not cover three host routes")
for row in production:
    for field in ("r1_paths", "r2_paths"):
        for path in row[field].split(";"):
            require(Path(path).is_file(), f"missing raw input: {path}")

snake_text = "\n".join(path.read_text() for path in [ROOT / "Snakefile", *sorted((ROOT / "workflow" / "rules").glob("*.smk"))])
required_rules = {
    "all", "phase1_ready", "trim_reads", "kraken_bracken", "sylph_screen",
    "phyloflash_screen", "index_host_reference", "host_handling", "validate_sample",
    "phase1_smoke", "screen_cohort",
}
for rule in required_rules:
    require(f"rule {rule}:" in snake_text, f"required rule missing: {rule}")

for environment in sorted((ROOT / "envs").glob("*.yaml")):
    text = environment.read_text()
    require("=" in text, f"environment lacks pinned dependencies: {environment.name}")
require(len(list((ROOT / "envs").glob("*.yaml"))) == 6, "expected six pinned environment definitions")

for path in [
    ROOT / "profiles" / "slurm" / "config.yaml",
    ROOT / "profiles" / "slurm" / "submit.py",
    ROOT / "profiles" / "slurm" / "status.py",
    ROOT / "scripts" / "workflow_controller.sbatch",
    ROOT / "scripts" / "submit_workflow.sh",
]:
    require(path.is_file(), f"missing cluster artifact: {path.relative_to(ROOT)}")

ready = ROOT / "data" / "results" / "stages" / "phase1.ready"
require(ready.is_file(), "production phase1.ready sentinel is missing")

smoke_root = ROOT / "tests" / "work"
smoke_stage = smoke_root / "stages" / "phase1_smoke.done"
validations = sorted((smoke_root / "validation" / "samples").glob("*.json"))
if args.require_smoke:
    require(smoke_stage.is_file(), "smoke stage sentinel is missing")
    require(len(validations) == 3, f"expected 3 smoke validations, found {len(validations)}")
    for path in validations:
        try:
            require(json.loads(path.read_text()).get("status") == "PASS", f"smoke validation failed: {path.name}")
        except Exception as exc:
            errors.append(f"cannot parse smoke validation {path}: {exc}")

metrics = {
    "status": "PASS" if not errors else "FAIL",
    "production_samples": len(production),
    "fixture_samples": len(fixture),
    "fixture_studies": len({row["study"] for row in fixture}),
    "fixture_host_routes": len({row["host_route"] for row in fixture}),
    "pinned_environments": len(list((ROOT / "envs").glob("*.yaml"))),
    "smoke_required": "yes" if args.require_smoke else "no",
    "smoke_validations": len(validations),
    "errors": len(errors),
}
REPORT.parent.mkdir(parents=True, exist_ok=True)
REPORT.write_text("".join(f"{key}\t{value}\n" for key, value in metrics.items()) +
                  "".join(f"ERROR\t{error}\n" for error in errors))
print(REPORT.read_text(), end="")
if errors:
    raise SystemExit(1)
