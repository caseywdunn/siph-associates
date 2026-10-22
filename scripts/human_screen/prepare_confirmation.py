"""Deduplicate a bounded, deterministic set of candidate pairs for broad nt search."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

from Bio import SeqIO


def prepare(root, samples, limit, output):
    output.mkdir(parents=True, exist_ok=True)
    selected, sequences = [], {}
    for sample in samples:
        directory = root / "samples" / sample
        records = {
            record.id: str(record.seq)
            for record in SeqIO.parse(directory / "candidates.fasta", "fasta")
        }
        with (directory / "pairs.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        # Select unique sequence pairs first; hash order is repeatable and unrelated
        # to alignment position. This is a candidate-validation sample, not abundance.
        unique = {}
        for row in rows:
            unique.setdefault(row["sequence_sha256"], row)
        for digest, row in sorted(unique.items())[:limit]:
            for mate in (1, 2):
                original = row["pair_id"] + f"/{mate}"
                sequence = records[original]
                query_id = "seq_" + hashlib.sha256(sequence.encode()).hexdigest()
                sequences[query_id] = sequence
                selected.append(
                    {
                        "sample_id": sample,
                        "pair_id": row["pair_id"],
                        "mate": mate,
                        "sequence_pair_sha256": digest,
                        "query_id": query_id,
                    }
                )
    with (output / "queries.fasta").open("w") as handle:
        for name, sequence in sorted(sequences.items()):
            handle.write(f">{name}\n{sequence}\n")
    with (output / "query_sources.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "pair_id",
                "mate",
                "sequence_pair_sha256",
                "query_id",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(selected)
    (output / "selection.json").write_text(
        json.dumps(
            {
                "samples": samples,
                "max_distinct_pairs_per_library": limit,
                "selection": "ascending_sha256_of_distinct_sequence_pairs",
                "unique_read_queries": len(sequences),
                "selected_pair_occurrences": len(selected) // 2,
            },
            indent=2,
        )
        + "\n"
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--scope", choices=("pilot", "cohort"), required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-pairs", type=int, default=500)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    if args.scope == "pilot":
        samples = config["pilot_samples"]
    else:
        with open(config["manifest"]) as handle:
            samples = [
                row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")
            ]
    prepare(Path(config["run"]), samples, args.max_pairs, args.output)


if __name__ == "__main__":
    main()
