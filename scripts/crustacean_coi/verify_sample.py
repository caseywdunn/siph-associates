#!/usr/bin/env python3
"""Apply the tested origin rules to one library's preserved COI candidates."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from verify import (
    decide_pair,
    hits,
    near_best,
    nucleotide_evidence,
    protein_evidence,
    table,
    write_table,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    curated = Path(config["run"]) / "curated"
    metadata = {r["reference_id"]: r for r in table(curated / "protein_metadata.tsv")}
    nucleotide_meta = {
        r["reference_id"]: r for r in table(curated / "reference_metadata.tsv")
    }
    protein = defaultdict(list)
    for code in [1, 5]:
        for query, records in hits(args.output / f"code{code}.tsv").items():
            protein[query].extend(records)
    nucleotide = defaultdict(list)
    with (args.directory / "competitive.tsv").open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            nucleotide[row[0]].append(
                dict(
                    reference_id=row[1],
                    identity=float(row[2]),
                    length=int(row[3]),
                    bitscore=float(row[8]),
                )
            )
    rows = []
    for row in table(args.directory / "candidate_pairs.tsv"):
        queries = [f"{row['candidate_id']}_{mate}" for mate in [1, 2]]
        aa = [protein_evidence(protein[q], metadata, config) for q in queries]
        nt = [
            nucleotide_evidence(nucleotide[q], nucleotide_meta, config) for q in queries
        ]
        refs = sorted(
            {
                h["reference_id"]
                for q in queries
                for h in near_best(protein[q], config["read_near_best_fraction"])
            }
        )
        rows.append(
            dict(
                **row,
                classification=decide_pair([a[0] for a in aa], nt),
                protein_R1=aa[0][0],
                protein_R2=aa[1][0],
                nucleotide_R1=nt[0],
                nucleotide_R2=nt[1],
                near_best_proteins=";".join(refs),
            )
        )
    if rows:
        write_table(args.output / "classified_pairs.tsv", rows)
    else:
        (args.output / "classified_pairs.tsv").write_text(
            "candidate_id\tread_pair_id\tbbduk\tdiamond\tclassification\tprotein_R1\tprotein_R2\tnucleotide_R1\tnucleotide_R2\tnear_best_proteins\n"
        )
    summary = json.loads((args.directory / "recruitment.json").read_text())
    summary.update(
        verification_config=str(args.config),
        classification_counts=Counter(r["classification"] for r in rows),
        interpretation="All recruited pairs retained; supported COI origin is not a species or dietary assignment.",
    )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
