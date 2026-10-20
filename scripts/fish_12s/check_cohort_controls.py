#!/usr/bin/env python3
"""Check paired recruitment and competitive classification of known controls."""

import argparse
import csv
import json
from pathlib import Path

from Bio import SeqIO

from classify_candidates import evidence, pair_category, read_hits, reference_metadata
from pilot import digest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--reference-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path)
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
    hits_path = out / "alignment/hits.tsv"
    hits = read_hits(hits_path, reference_metadata(args.reference_metadata))
    if set(hits) - query_ids:
        raise ValueError("Control alignment contains unknown query identifiers")
    classified = {
        r["read_pair_id"]: pair_category(
            {
                mate: evidence(hits.get(f"{r['candidate_id']}_{mate}", []))
                for mate in [1, 2]
            }
        )
        for r in rows
    }
    expected_fish = {
        "motif_R1",
        "motif_reverse_R2",
        "fish_away_from_motif",
        "cartilaginous_motif_variant",
    }
    if {
        rid
        for rid, category in classified.items()
        if category == "long_close_fish_match"
    } != expected_fish:
        raise ValueError(f"Competitive control classification failed: {classified}")
    if classified.get("frog") != "other_or_no_hit":
        raise ValueError("The frog control was not rejected as non-fish")
    frog = next(r for r in rows if r["read_pair_id"] == "frog")
    frog_hit = evidence(hits.get(f"{frog['candidate_id']}_1", []))["top"]
    if (
        frog_hit.get("fish", True)
        or frog_hit.get("identity", 0) < 99
        or frog_hit.get("length", 0) < 150
    ):
        raise ValueError("The frog control lacks its expected close non-fish alignment")
    marker = args.output or args.root / "CONTROLS_PASSED"
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(
        json.dumps(
            dict(
                passed=True,
                input_pairs=9,
                candidate_pairs=5,
                exported_reads=10,
                expected_pair_ids=sorted(expected),
                classifications=classified,
                competitive_classification_checked=True,
                input_sha256={
                    str(p): digest(p)
                    for p in [
                        hits_path,
                        args.reference_metadata,
                        args.truth,
                        out / "candidates.fasta",
                        out / "candidate_pairs.tsv",
                        out / "run.json",
                    ]
                },
                source_sha256={
                    p.name: digest(p)
                    for p in [
                        Path(__file__),
                        Path(__file__).with_name("classify_candidates.py"),
                    ]
                },
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
