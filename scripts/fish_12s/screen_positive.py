#!/usr/bin/env python3
"""Recover 12S candidates from complete lanes selected as positive controls."""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from pilot import MOTIF, digest, pair_ids, read_id, records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    with args.manifest.open() as handle:
        row = list(csv.DictReader(handle, delimiter="\t"))[args.index]
    out = args.output_root / row["analysis_sample_id"] / row["read_pair_id"]
    out.mkdir(parents=True, exist_ok=False)
    inputs = [Path(row[f"r{mate}_path"]) for mate in (1, 2)]
    for mate, path in enumerate(inputs, 1):
        if path.stat().st_size != int(row[f"r{mate}_bytes"]):
            raise ValueError(f"Input size differs from inventory: {path}")
    provenance = {
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "status": "running",
        "command": sys.argv,
        "python": sys.version,
        "input": row,
        "raw_input_stat": {
            str(p): {"size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
            for p in inputs
        },
        "manifest_sha256": digest(args.manifest),
        "reference": str(args.reference),
        "reference_sha256": digest(args.reference),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "seqkit_version": subprocess.check_output(
            ["seqkit", "version"], text=True
        ).strip(),
        "source_sha256": {},
        "commands": [],
    }
    for source in [Path(__file__), Path(__file__).with_name("pilot.py")]:
        shutil.copy2(source, out / source.name)
        provenance["source_sha256"][source.name] = digest(source)

    def save() -> None:
        (out / "run.json").write_text(json.dumps(provenance, indent=2) + "\n")

    def run(name: str, command: list[str], stdout: Path | None = None) -> None:
        start = time.monotonic()
        with (
            (out / f"{name}.log").open("wb") as log,
            (stdout or out / f"{name}.stdout").open("wb") as output,
        ):
            result = subprocess.run(command, stdout=output, stderr=log, check=False)
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

    def bbduk(k: int, paths: list[Path], label: str) -> set[str]:
        run(
            label,
            [
                "bbduk.sh",
                "-Xmx4g",
                f"threads={args.threads}",
                f"k={k}",
                "hdist=0",
                "maskmiddle=f",
                "rcomp=t",
                "rieb=t",
                "qtrim=f",
                "minlen=0",
                "ordered=t",
                "overwrite=f",
                "pigz=f",
                "unpigz=f",
                f"in={paths[0]}",
                f"in2={paths[1]}",
                f"ref={args.reference}",
                f"outm={out / (label + '_R1.fastq')}",
                f"outm2={out / (label + '_R2.fastq')}",
            ],
        )
        return pair_ids(out / f"{label}_R1.fastq", out / f"{label}_R2.fastq")

    save()
    try:
        k25 = bbduk(25, inputs, "bbduk25")
        # Every exact 31-mer contains exact 25-mers. This nested recheck has
        # identical candidate membership to a full k=31 scan, but is not a
        # measurement of full-input k=31 runtime.
        k31 = bbduk(31, [out / f"bbduk25_R{m}.fastq" for m in (1, 2)], "bbduk31")
        if not k31 <= k25:
            raise ValueError(
                "Exact 31-mer candidates must be a subset of 25-mer candidates"
            )
        motif_ids, matched_reads = set(), set()
        for mate, path in enumerate(inputs, 1):
            run(
                f"regex_R{mate}",
                [
                    "seqkit",
                    "grep",
                    "-j",
                    str(args.threads),
                    "-s",
                    "-r",
                    "-P",
                    "-p",
                    MOTIF,
                    str(path),
                ],
                stdout=out / f"regex_R{mate}.fastq",
            )
            for name, seq, _ in records(out / f"regex_R{mate}.fastq"):
                if re.search(MOTIF, seq) is None:
                    raise ValueError(f"Regex output failed sequence-only check: {name}")
                motif_ids.add(read_id(name))
                matched_reads.add((read_id(name), mate))

        # Preserve raw IDs and qualities in the FASTQ outputs. Query IDs in the
        # small alignment FASTA map back to their source in candidate_reads.tsv.
        seqs = {}
        for label in ("bbduk25", "regex"):
            for mate in (1, 2):
                for name, seq, qual in records(out / f"{label}_R{mate}.fastq"):
                    key = (read_id(name), mate)
                    if key in seqs and seqs[key] != (name, seq, qual):
                        raise ValueError(f"Conflicting candidate records: {key}")
                    seqs[key] = (name, seq, qual)
        with (
            (out / "candidates.fasta").open("w") as fasta,
            (out / "candidate_reads.tsv").open("w") as table,
        ):
            writer = csv.writer(table, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "query_id",
                    "read_pair_id",
                    "mate",
                    "raw_header",
                    "regex",
                    "bbduk25_pair",
                    "bbduk31_pair",
                ]
            )
            for index, ((rid, mate), (name, seq, _)) in enumerate(
                sorted(seqs.items()), 1
            ):
                qid = f"q{index:07d}"
                fasta.write(f">{qid}\n{seq}\n")
                writer.writerow(
                    [
                        qid,
                        rid,
                        mate,
                        name,
                        int((rid, mate) in matched_reads),
                        int(rid in k25),
                        int(rid in k31),
                    ]
                )
        match = re.search(r"Input:\s+(\d+) reads", (out / "bbduk25.log").read_text())
        if match is None or int(match[1]) % 2:
            raise ValueError("Cannot recover the paired input count from BBDuk's log")
        summary = {
            "analysis_sample_id": row["analysis_sample_id"],
            "read_pair_id": row["read_pair_id"],
            "pairs_screened": int(match[1]) // 2,
            "regex_reads": len(matched_reads),
            "regex_pairs": len(motif_ids),
            "bbduk25_pairs": len(k25),
            "bbduk31_pairs": len(k31),
            "regex_and_bbduk25_pairs": len(motif_ids & k25),
            "regex_and_bbduk31_pairs": len(motif_ids & k31),
            "regex_only_pairs": len(motif_ids - k25),
            "bbduk25_not_regex_pairs": len(k25 - motif_ids),
            "candidate_reads": len(seqs),
        }
        (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        provenance["status"] = "complete"
        provenance["outputs"] = {
            p.name: digest(p)
            for p in out.iterdir()
            if p.is_file() and p.name != "run.json"
        }
        print(json.dumps(summary))
    except Exception as exc:
        provenance["status"] = "failed"
        provenance["error"] = repr(exc)
        raise
    finally:
        provenance["finished_utc"] = datetime.now(timezone.utc).isoformat()
        save()


if __name__ == "__main__":
    main()
