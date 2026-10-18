#!/usr/bin/env python3
"""Check complete paired export against the pilot's known-answer fish controls."""

import argparse
import csv
import json
from pathlib import Path

from Bio import SeqIO


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    args = parser.parse_args()
    out = args.root / "union/synthetic/host_depleted"
    expected = set(json.loads(args.truth.read_text())["expected_pair_ids"])
    with (out / "candidate_pairs.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    observed = {r["read_pair_id"] for r in rows}
    if observed != expected or len(rows) != 5:
        raise ValueError("Unexpected fish control recruitment")
    query_ids = {r.id for r in SeqIO.parse(out / "candidates.fasta", "fasta")}
    if query_ids != {f"{r['candidate_id']}_{mate}" for r in rows for mate in [1, 2]}:
        raise ValueError("Each control must export both original mates")
    summary = json.loads((out / "run.json").read_text())
    if summary["input_pairs"] != 9 or summary["method_pairs"] != dict(
        regex=3, bbduk25=5, sortmerna=5, bowtie2=5
    ):
        raise ValueError("Control accounting differs from the validated pilot")
    (args.root / "CONTROLS_PASSED").write_text(
        json.dumps(
            dict(
                passed=True,
                input_pairs=9,
                candidate_pairs=5,
                exported_reads=10,
                expected_pair_ids=sorted(expected),
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
