"""Summarize human-compatible candidates against broad nucleotide alternatives."""

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

FIELDS = "qseqid saccver pident length qlen qstart qend evalue bitscore staxids sscinames".split()
UNKNOWN_TAXA = {"", "N/A", "0"}


def classify_query(rows, subject_count):
    """Keep truncated or taxonomically incomplete candidate searches unresolved."""
    qualified = [
        r
        for r in rows
        if float(r["pident"]) >= 98
        and (int(r["qend"]) - int(r["qstart"]) + 1) / int(r["qlen"]) >= 0.95
    ]
    best = max((float(r["bitscore"]) for r in qualified), default=0)
    near = [r for r in qualified if float(r["bitscore"]) >= 0.95 * best]
    taxa = {taxon for row in near for taxon in row["staxids"].split(";")}
    names = {name for row in near for name in row["sscinames"].split(";")}
    human = "9606" in taxa
    other = bool(taxa - UNKNOWN_TAXA - {"9606"})
    unknown = bool(taxa & UNKNOWN_TAXA)
    capped = subject_count >= 100
    if capped or unknown:
        classification = "unresolved_truncation_or_unknown_taxonomy"
    elif not near:
        classification = "no_qualifying_hit"
    elif human and other:
        classification = "human_with_alternatives"
    elif human:
        classification = "human_only_near_best"
    else:
        classification = "human_not_near_best"
    return {
        "best_bitscore": best,
        "near_best_scientific_names": ";".join(sorted(names)),
        "near_best_taxids": ";".join(sorted(taxa)),
        "human_among_near_best": human,
        "nonhuman_among_near_best": other,
        "unknown_near_best_taxonomy": unknown,
        "hit_limit_reached": capped,
        "query_classification": classification,
    }


def classify_pair(mates):
    """Require both mates to support a pair-level taxonomic conclusion."""
    if len(mates) != 2:
        raise ValueError("Expected two candidate mates")
    labels = [mate["query_classification"] for mate in mates]
    if "unresolved_truncation_or_unknown_taxonomy" in labels:
        return "unresolved_truncation_or_unknown_taxonomy"
    if "no_qualifying_hit" in labels:
        return "incomplete_confirmation"
    if all(label == "human_only_near_best" for label in labels):
        return "both_mates_human_only_near_best"
    if all(mate["human_among_near_best"] for mate in mates):
        return "human_with_alternatives"
    if any(mate["human_among_near_best"] for mate in mates):
        return "conflicting_mates"
    return "human_not_supported_by_nt"


def scientific_names(path, wanted):
    """Resolve existing taxids without rerunning the sequence search."""
    names = {}
    with path.open() as handle:
        for line in handle:
            fields = [field.strip() for field in line.split("|")]
            if fields[0] in wanted and fields[3] == "scientific name":
                names[fields[0]] = fields[1]
                if len(names) == len(wanted):
                    break
    return names


def review(hits_path, sources_path, output, names_path=None):
    hits, subjects = defaultdict(list), defaultdict(set)
    with hits_path.open() as handle:
        for row in csv.DictReader(handle, fieldnames=FIELDS, delimiter="\t"):
            hits[row["qseqid"]].append(row)
            subjects[row["qseqid"]].add(row["saccver"])
    if names_path:
        wanted = {
            taxon
            for rows in hits.values()
            for row in rows
            for taxon in row["staxids"].split(";")
        } - UNKNOWN_TAXA
        names = scientific_names(names_path, wanted)
        for rows in hits.values():
            for row in rows:
                if row["sscinames"] in ("", "N/A"):
                    row["sscinames"] = ";".join(
                        names.get(taxon, "taxid:" + taxon)
                        for taxon in row["staxids"].split(";")
                    )
    with sources_path.open() as handle:
        sources = list(csv.DictReader(handle, delimiter="\t"))
    source_fields = ["sample_id", "pair_id", "mate", "sequence_pair_sha256", "query_id"]
    review_fields = list(classify_query([], 0))
    pair_mates, query_counts = defaultdict(list), Counter()
    with output.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=source_fields + review_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for source in sources:
            query = source["query_id"]
            assessed = classify_query(hits[query], len(subjects[query]))
            query_counts[assessed["query_classification"]] += 1
            writer.writerow({**source, **assessed})
            pair_mates[(source["sample_id"], source["pair_id"])].append(assessed)
    pair_counts, sample_counts = Counter(), defaultdict(Counter)
    with output.with_name("pair_review.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "pair_id",
                "pair_classification",
                "mate1_classification",
                "mate2_classification",
            ]
        )
        for (sample, pair), mates in sorted(pair_mates.items()):
            label = classify_pair(mates)
            pair_counts[label] += 1
            sample_counts[sample][label] += 1
            writer.writerow(
                [sample, pair, label, *(m["query_classification"] for m in mates)]
            )
    output.with_suffix(".json").write_text(
        json.dumps(
            {
                "selected_pair_occurrences": len(pair_mates),
                "pair_classification_counts": dict(pair_counts),
                "mate_classification_counts": dict(query_counts),
                "per_sample_pair_classifications": {
                    key: dict(value) for key, value in sample_counts.items()
                },
            },
            indent=2,
        )
        + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", required=True, type=Path)
    parser.add_argument("--sources", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--names", type=Path)
    args = parser.parse_args()
    review(args.hits, args.sources, args.output, args.names)


if __name__ == "__main__":
    main()
