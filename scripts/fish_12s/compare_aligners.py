#!/usr/bin/env python3
"""Compare direct 12S alignment filters on the same complete paired inputs."""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from pilot import digest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--bowtie-index", type=Path, required=True)
    parser.add_argument("--sortmerna", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    with args.manifest.open() as handle:
        row = list(csv.DictReader(handle, delimiter="\t"))[args.index]
    out = args.output_root / row["analysis_sample_id"] / row["read_pair_id"]
    out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__)
    shutil.copy2(source, out / source.name)
    shutil.copy2(source.with_name("pilot.py"), out / "pilot.py")
    provenance = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
        "command": sys.argv,
        "python": sys.version,
        "input": row,
        "manifest_sha256": digest(args.manifest),
        "reference_sha256": digest(args.reference),
        "source_sha256": {p.name: digest(p) for p in out.glob("*.py")},
        "bowtie_index_sha256": {
            p.name: digest(p)
            for p in args.bowtie_index.parent.glob(args.bowtie_index.name + ".*.bt2")
        },
        "input_stat": {
            row[f"r{m}_path"]: {
                "size": Path(row[f"r{m}_path"]).stat().st_size,
                "mtime_ns": Path(row[f"r{m}_path"]).stat().st_mtime_ns,
            }
            for m in (1, 2)
        },
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "commands": [],
    }

    def save() -> None:
        (out / "run.json").write_text(json.dumps(provenance, indent=2) + "\n")

    def run(name: str, command: list[str]) -> None:
        start = time.monotonic()
        with (out / f"{name}.log").open("wb") as log:
            result = subprocess.run(
                command, stdout=log, stderr=subprocess.STDOUT, check=False
            )
        provenance["commands"].append(
            {
                "name": name,
                "command": command,
                "seconds": time.monotonic() - start,
                "returncode": result.returncode,
            }
        )
        save()
        result.check_returncode()

    save()
    try:
        run("bowtie2_version", ["bowtie2", "--version"])
        run("sortmerna_version", [str(args.sortmerna), "--version"])
        run(
            "bowtie2",
            [
                "bowtie2",
                "--very-sensitive-local",
                "--seed",
                "20261017",
                "--no-unal",
                "-p",
                str(args.threads),
                "-x",
                str(args.bowtie_index),
                "-1",
                row["r1_path"],
                "-2",
                row["r2_path"],
                "-S",
                str(out / "bowtie2.sam"),
            ],
        )
        run(
            "sortmerna",
            [
                str(args.sortmerna),
                "--ref",
                str(args.reference),
                "--reads",
                row["r1_path"],
                "--reads",
                row["r2_path"],
                "--threads",
                str(args.threads),
                "--workdir",
                str((out / "sortmerna_work").resolve()),
                "--aligned",
                str((out / "sortmerna").resolve()),
                "--fastx",
                "--paired_in",
                "--zip-out",
                "false",
                "--blast",
                "1",
            ],
        )
        provenance["status"] = "complete"
        provenance["outputs"] = {
            p.name: digest(p)
            for p in out.iterdir()
            if p.is_file() and p.name != "run.json"
        }
    except Exception as exc:
        provenance["status"] = "failed"
        provenance["error"] = repr(exc)
        raise
    finally:
        provenance["finished_utc"] = datetime.now(timezone.utc).isoformat()
        save()


if __name__ == "__main__":
    main()
