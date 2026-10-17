#!/usr/bin/env python3
"""Compare the original regex, record-aware matching, and whole-12S k-mers."""

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from itertools import zip_longest
from pathlib import Path

from Bio.SeqIO.QualityIO import FastqGeneralIterator

MOTIF = "ACAAAC..GGATT..ATACCCC..TATGC|GCATA..GGGGTAT..AATCC..GTTTGT"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        yield from FastqGeneralIterator(handle)


def read_id(name: str) -> str:
    return re.sub(r"/[12]$", "", name.split()[0])


def pair_ids(first: Path, second: Path) -> set[str]:
    ids = set()
    for a, b in zip_longest(records(first), records(second)):
        if a is None or b is None or read_id(a[0]) != read_id(b[0]):
            raise ValueError(f"Unsynchronized mates: {first}, {second}")
        rid = read_id(a[0])
        if rid in ids:
            raise ValueError(f"Duplicate read pair: {rid}")
        ids.add(rid)
    return ids


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--pairs", type=int, default=1000000)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--bbduk", required=True)
    args = parser.parse_args()
    with args.manifest.open() as handle:
        row = list(csv.DictReader(handle, delimiter="\t"))[args.index]
    out = args.output_root / row["analysis_sample_id"]
    out.mkdir(parents=True, exist_ok=False)
    started = datetime.now(timezone.utc).isoformat()
    commands = []

    def run(
        name: str,
        command: list[str],
        stdout: Path | None = None,
        env: dict | None = None,
        empty_ok: bool = False,
    ) -> None:
        begin = time.monotonic()
        with (
            (out / f"{name}.log").open("wb") as log,
            (stdout or out / f"{name}.stdout").open("wb") as output,
        ):
            result = subprocess.run(
                command, stdout=output, stderr=log, env=env, check=False
            )
        entry = {
            "name": name,
            "command": command,
            "seconds": time.monotonic() - begin,
            "returncode": result.returncode,
        }
        if env is not None:
            entry["LC_ALL"] = env.get("LC_ALL", "")
        commands.append(entry)
        (out / "commands.json").write_text(json.dumps(commands, indent=2) + "\n")
        if result.returncode != 0 and not (empty_ok and result.returncode == 1):
            raise RuntimeError(f"{name} failed: see {out / (name + '.log')}")

    # The same bounded prefix is used by every method, with no random sampling.
    paths = [out / "input_R1.fastq.gz", out / "input_R2.fastq.gz"]
    for mate, path in enumerate(paths, 1):
        run(
            f"subset_R{mate}",
            [
                "seqkit",
                "head",
                "-j",
                "1",
                "-n",
                str(args.pairs),
                "-o",
                str(path),
                row[f"r{mate}_path"],
            ],
        )
    original_ids = pair_ids(*paths)
    if len(original_ids) != args.pairs:
        raise ValueError("Source lane contains fewer pairs than requested")
    for mate, path in enumerate(paths, 1):
        run(
            f"zgrep_R{mate}",
            ["zgrep", "-E", MOTIF, str(path)],
            stdout=out / f"zgrep_R{mate}.txt",
            empty_ok=True,
        )
        run(
            f"zgrep_C_R{mate}",
            ["zgrep", "-E", MOTIF, str(path)],
            stdout=out / f"zgrep_C_R{mate}.txt",
            env={**os.environ, "LC_ALL": "C"},
            empty_ok=True,
        )
        run(
            f"seqkit_R{mate}",
            ["seqkit", "grep", "-j", "1", "-s", "-r", "-P", "-p", MOTIF, str(path)],
            stdout=out / f"motif_R{mate}.fastq",
        )
    motif_ids = set()
    motif_sequences = []
    for mate in [1, 2]:
        seqs = [seq for _, seq, _ in records(out / f"motif_R{mate}.fastq")]
        raw_lines = (out / f"zgrep_R{mate}.txt").read_text().splitlines()
        c_lines = (out / f"zgrep_C_R{mate}.txt").read_text().splitlines()
        if raw_lines != c_lines or raw_lines != seqs:
            raise ValueError(
                "grep and sequence-only matches differ; inspect quality/header matches"
            )
        motif_sequences.extend(seqs)
        motif_ids.update(
            read_id(name) for name, _, _ in records(out / f"motif_R{mate}.fastq")
        )
    # Explicit exact 31-mers; no quality/length trimming or middle-base masking.
    run(
        "bbduk",
        [
            args.bbduk,
            "-Xmx4g",
            "threads=4",
            "k=31",
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
            f"outm={out / 'kmer_R1.fastq'}",
            f"outm2={out / 'kmer_R2.fastq'}",
        ],
    )
    kmer_ids = pair_ids(out / "kmer_R1.fastq", out / "kmer_R2.fastq")
    if not (motif_ids | kmer_ids) <= original_ids:
        raise ValueError("Filter returned IDs not in the input")
    summary = {
        "analysis_sample_id": row["analysis_sample_id"],
        "pairs_screened": len(original_ids),
        "motif_reads": len(motif_sequences),
        "motif_pairs": len(motif_ids),
        "whole_12s_kmer_pairs": len(kmer_ids),
        "both_pairs": len(motif_ids & kmer_ids),
        "motif_only_pairs": len(motif_ids - kmer_ids),
        "kmer_only_pairs": len(kmer_ids - motif_ids),
    }
    # Preserve both mates and qualities for the union of both candidate filters.
    with (
        (out / "candidate_R1.fastq").open("w") as f1,
        (out / "candidate_R2.fastq").open("w") as f2,
    ):
        for a, b in zip(records(paths[0]), records(paths[1])):
            if read_id(a[0]) in (motif_ids | kmer_ids):
                for handle, (name, seq, qual) in [(f1, a), (f2, b)]:
                    handle.write(f"@{name}\n{seq}\n+\n{qual}\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    provenance = {
        "started_utc": started,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv,
        "python": sys.version,
        "script_sha256": digest(Path(__file__)),
        "manifest_sha256": digest(args.manifest),
        "reference_sha256": digest(args.reference),
        "input": row,
        "raw_input_stat": {
            str(p): {"size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
            for p in [Path(row["r1_path"]), Path(row["r2_path"])]
        },
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "seqkit_version": subprocess.check_output(
            ["seqkit", "version"], text=True
        ).strip(),
        "grep_version": subprocess.check_output(
            ["grep", "--version"], text=True
        ).splitlines()[0],
        "outputs": {str(p): digest(p) for p in out.iterdir() if p.is_file()},
    }
    (out / "run.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
