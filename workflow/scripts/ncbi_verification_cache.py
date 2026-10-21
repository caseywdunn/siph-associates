"""Reuse dated NCBI verification only for exactly identical query sequences."""

import csv
import hashlib
from collections import defaultdict
from pathlib import Path


def fasta_sequences(path):
    records, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in records:
                raise ValueError(f"Duplicate FASTA identifier: {name}")
            records[name] = ""
        elif name:
            records[name] += line.strip().upper()
    return records


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def cached_verifications(queries, fasta, table):
    """Return all recorded hits for matching sequences, retaining original query dates/RIDs."""
    originals = fasta_sequences(fasta)
    groups = defaultdict(list)
    with open(table, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sequence_id"] not in originals:
                raise ValueError("Verification row is absent from its source FASTA")
            groups[row["sequence_id"]].append(row)
    by_sequence = {}
    for name, hits in sorted(
        groups.items(),
        key=lambda item: (max(r["query_date"] for r in item[1]), item[0]),
    ):
        by_sequence[originals[name]] = hits
    reused = []
    pending = []
    for name, sequence in queries.items():
        if sequence not in by_sequence:
            pending.append(name)
            continue
        for hit in by_sequence[sequence]:
            reused.append(
                {
                    **hit,
                    "sequence_id": name,
                    "sequence_sha256": sequence_sha256(sequence),
                    "verification_source": f"{table}#{hit['sequence_id']}",
                }
            )
    return reused, pending
