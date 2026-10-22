#!/usr/bin/env python3
"""Run one bounded community sensitivity and record its exact inputs."""

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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(config_path: Path) -> None:
    started = datetime.now(timezone.utc)
    settings = json.loads(config_path.read_text())
    inputs = [
        Path(settings[key])
        for key in ("metadata", "membership", "original_tests", "bacterial_grades")
    ]
    sources = [
        config_path,
        Path(settings["original_config"]),
        Path("workflow/scripts/analyze_coauthor_communities.R"),
        Path("workflow/scripts/analyze_cutibacterium_sensitivity.R"),
        Path("scripts/run_cutibacterium_sensitivity.py"),
        Path("workflow/cutibacterium_sensitivity.smk"),
        Path("envs/stats.yaml"),
    ]
    original_hashes = {str(path): sha256(path) for path in inputs + sources}
    command = [
        settings["rscript"],
        "workflow/scripts/analyze_cutibacterium_sensitivity.R",
        str(config_path),
    ]
    subprocess.run(command, check=True)
    if original_hashes != {str(path): sha256(path) for path in inputs + sources}:
        raise ValueError("Inputs or sources changed during the sensitivity calculation")
    out = Path(settings["output_dir"])
    for source in sources:
        destination = out / "provenance/source" / source
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    outputs = [
        out / name
        for name in (
            "community_tests.tsv",
            "community_membership.tsv",
            "summary.json",
            "values.tex",
            "session_info.txt",
        )
    ]
    result = json.loads((out / "summary.json").read_text())
    provenance = {
        "analysis_id": "community-cutibacterium-omission-sensitivity",
        "run_id": started.strftime("%Y%m%dT%H%M%SZ"),
        "started_at_utc": started.isoformat(),
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "command": shlex.join([sys.executable, *sys.argv]),
        "r_command": shlex.join(command),
        "python_version": platform.python_version(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "base_git_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True
        ).strip(),
        "inputs_sha256": {str(path): original_hashes[str(path)] for path in inputs},
        "source_sha256": {str(path): original_hashes[str(path)] for path in sources},
        "source_snapshot": str(out / "provenance/source"),
        "outputs_sha256": {str(path): sha256(path) for path in outputs},
        "checks": {
            "unchanged_inputs_and_sources": True,
            "baseline_reconstruction_matches": result[
                "baseline_reconstruction_matches"
            ],
        },
        "parameters": settings,
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/cutibacterium_sensitivity.json")
    )
    run(parser.parse_args().config)
