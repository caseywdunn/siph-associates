#!/usr/bin/env python3
"""Compare protein-level competitive assignment on the held-out COI controls."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

from summarize import classify, combine


def prepare(root: Path, source: Path) -> None:
    metadata = {
        r["reference_id"]: r
        for r in csv.DictReader(
            (root / "panel/reference_metadata.tsv").open(), delimiter="\t"
        )
    }
    output = root / "protein_competition"
    output.mkdir(parents=True, exist_ok=True)
    used = set()
    with (output / "reference.faa").open("w") as handle:
        for path in [root / "panel/target.faa", source / "annotated.faa"]:
            for record in SeqIO.parse(path, "fasta"):
                if record.id in metadata and record.id not in used:
                    SeqIO.write(record, handle, "fasta")
                    used.add(record.id)
    (output / "reference_summary.json").write_text(
        json.dumps(
            dict(
                proteins=len(used),
                target=sum(metadata[r]["target"] == "1" for r in used),
                non_target=sum(metadata[r]["target"] == "0" for r in used),
                scope="Withheld panel target proteins plus retained annotated RefSeq COI competitors; unannotated cnidarian barcodes and bacterial COX homologs are not represented.",
            ),
            indent=2,
        )
        + "\n"
    )


def report(root: Path, all_near_best: bool = False) -> None:
    metadata = {
        r["reference_id"]: r
        for r in csv.DictReader(
            (root / "panel/reference_metadata.tsv").open(), delimiter="\t"
        )
    }
    pairs = {
        r["read_pair_id"]: r
        for r in csv.DictReader(
            (root / "heldout/candidate_pairs.tsv").open(), delimiter="\t"
        )
    }
    truth = json.loads((root / "heldout/truth.json").read_text())
    hits = defaultdict(list)
    name = "all_near_best_hits.tsv" if all_near_best else "hits.tsv"
    with (root / "protein_competition" / name).open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            hits[row[0]].append(
                dict(
                    reference_id=row[1],
                    identity=float(row[2]),
                    length=int(row[3]) * 3,
                    bitscore=float(row[4]),
                )
            )
    summaries, details = {}, []
    for threshold in [60, 70, 80]:
        counts = defaultdict(Counter)
        for row in truth:
            pair = pairs.get(row["read_pair_id"])
            if pair:
                # Reuse the original score ranking, tie, truncation and mate rules.
                mates = []
                refs = set()
                for mate in [1, 2]:
                    query = f"{pair['candidate_id']}_{mate}"
                    adjusted = [
                        dict(h, identity=80 if h["identity"] >= threshold else 79)
                        for h in hits[query]
                    ]
                    label, _, references = classify(
                        adjusted, metadata, max_hits=None if all_near_best else 200
                    )
                    mates.append(label)
                    refs.update(references)
                label = combine(mates)
            else:
                label, refs = "not_recruited", set()
            counts[row["group"]]["pairs"] += 1
            counts[row["group"]][label] += 1
            details.append(
                dict(
                    read_pair_id=row["read_pair_id"],
                    expected=row["expected"],
                    group=row["group"],
                    protein_identity_threshold=threshold,
                    observed=label,
                    near_best_references=";".join(sorted(refs)),
                )
            )
        summaries[threshold] = counts
    output = root / "protein_competition"
    prefix = "all_near_best_" if all_near_best else ""
    (output / f"{prefix}summary.json").write_text(
        json.dumps(summaries, indent=2) + "\n"
    )
    with (output / f"{prefix}classified_controls.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle, list(details[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(details)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "report"])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--all-near-best", action="store_true")
    args = parser.parse_args()
    if args.stage == "prepare":
        prepare(args.root, args.source)
    else:
        report(args.root, args.all_near_best)


if __name__ == "__main__":
    main()
