#!/usr/bin/env python3
"""Checksum the Phase-1 workflow implementation and locked configuration."""
from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUTPUT = ROOT / "data" / "metadata" / "phase1.freeze.sha256"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


paths = [
    ROOT / "Snakefile",
    ROOT / "config" / "config.yaml",
    ROOT / "config" / "samples.tsv",
    ROOT / "tests" / "config.fixture.yaml",
    ROOT / "tests" / "fixtures" / "samples.tsv",
    ROOT / "scripts" / "build_workflow_samples.py",
    ROOT / "scripts" / "build_fixture_samples.py",
    ROOT / "scripts" / "validate_phase1.py",
    ROOT / "scripts" / "freeze_phase1.py",
    ROOT / "scripts" / "submit_workflow.sh",
    ROOT / "scripts" / "workflow_controller.sbatch",
] + sorted((ROOT / "envs").glob("*.yaml")) \
  + sorted((ROOT / "profiles" / "slurm").glob("*")) \
  + sorted((ROOT / "workflow" / "rules").glob("*.smk")) \
  + sorted((ROOT / "workflow" / "scripts").glob("*.py"))

OUTPUT.write_text("".join(
    f"{sha256(path)}  {path.relative_to(ROOT)}\n" for path in paths if path.is_file()
))
print(f"froze {len(paths)} Phase-1 files -> {OUTPUT}")
