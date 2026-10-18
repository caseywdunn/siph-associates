#!/usr/bin/env python3
"""Conservatively summarize competitive COI alignments; no species/prey calls."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


def classify(hits, metadata, max_hits=200):
    """Rank all competitors before applying the minimum evidence thresholds."""
    if not hits:
        return "unresolved", [], []
    best = max(h["bitscore"] for h in hits)
    near = [h for h in hits if h["bitscore"] >= 0.99 * best]
    ids = sorted({h["reference_id"] for h in near})
    target = {metadata[r]["target"] == "1" for r in ids}
    groups = sorted({metadata[r]["group"] for r in ids})
    # A full result list cannot establish that other near-best competitors were absent.
    truncated_near_best = (
        max_hits is not None
        and len(hits) >= max_hits
        and min(h["bitscore"] for h in hits) >= 0.99 * best
    )
    if truncated_near_best or len(target) > 1:
        return "ambiguous", groups, ids
    # A short, weaker match must not replace a stronger divergent competitor.
    if any(h["length"] < 90 or h["identity"] < 80 for h in near):
        return "unresolved", groups, ids
    return ("crustacean" if True in target else "non_target"), groups, ids


def combine(labels):
    informative = set(labels) - {"unresolved"}
    if "ambiguous" in informative or {"crustacean", "non_target"} <= informative:
        return "ambiguous"
    return next(iter(informative), "unresolved")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, help="Separate revised reports from preserved alignments"
    )
    args = parser.parse_args()
    root = args.directory
    output = args.output or root
    output.mkdir(parents=True, exist_ok=True)
    with args.metadata.open() as handle:
        metadata = {
            r["reference_id"]: r for r in csv.DictReader(handle, delimiter="\t")
        }
    hits = defaultdict(list)
    with (root / "competitive.tsv").open() as handle:
        for line in handle:
            query, ref, identity, length, qstart, qend, rstart, rend, score, evalue = (
                line.split()
            )
            hits[query].append(
                {
                    "reference_id": ref,
                    "identity": float(identity),
                    "length": int(length),
                    "bitscore": float(score),
                }
            )
    with (root / "candidate_pairs.tsv").open() as handle:
        candidates = list(csv.DictReader(handle, delimiter="\t"))
    counts = defaultdict(Counter)
    fields = [
        "candidate_id",
        "read_pair_id",
        "bbduk",
        "diamond",
        "classification",
        "groups",
        "near_best_species",
        "near_best_references",
    ]
    with (output / "classified_pairs.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in candidates:
            mates = [
                classify(hits.get(f"{row['candidate_id']}_{m}", []), metadata)
                for m in (1, 2)
            ]
            label = combine([m[0] for m in mates])
            references = sorted({r for m in mates for r in m[2]})
            species = sorted({metadata[r]["species"] for r in references})
            row.update(
                classification=label,
                groups=";".join(sorted({g for m in mates for g in m[1]})),
                near_best_species=";".join(species),
                near_best_references=";".join(references),
            )
            writer.writerow(row)
            counts["union"][label] += 1
            for method in ["bbduk", "diamond"]:
                if row[method] == "1":
                    counts[method][label] += 1
            if row["bbduk"] == row["diamond"] == "1":
                counts["both_methods"][label] += 1
    summary = json.loads((root / "recruitment.json").read_text())
    summary.update(
        classification_by_method=dict(counts),
        interpretation="Provisional COI similarity evidence, not validated species identification or proof of predation.",
        thresholds={
            "ranking": "All hits ranked before minimum length/identity filtering",
            "minimum_aligned_nt": 90,
            "minimum_identity_percent": 80,
            "near_best_bitscore_fraction": 0.99,
            "max_target_sequences": 200,
        },
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
