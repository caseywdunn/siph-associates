#!/usr/bin/env python3
"""Collect non-host assembled eukaryotic SSUs (phyloFlash and assembly 18S) for grading and verification."""
import argparse
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eukaryote_units import unit  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--classifications", required=True, nargs="+", help="<sample>.<source>.tsv LCA tables")
parser.add_argument("--fastas", required=True, nargs="+", help="matching query FASTA files, same order")
parser.add_argument("--table", required=True, type=Path)
parser.add_argument("--fasta", required=True, type=Path)
args = parser.parse_args()

rows, sequences = [], {}
for table, fasta in zip(args.classifications, args.fastas):
    source = Path(table).name.split(".")[-2]
    wanted = {}
    with open(table, newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            u = unit(r["lineage"])
            if u and u[0] in ("metazoan", "non_metazoan", "human"):
                key = f"{r['sample_id']}|{source}|{r['sequence_id']}"
                wanted[r["sequence_id"]] = key
                rows.append({"sample_id": r["sample_id"], "source": source, "sequence_id": key,
                             "lineage": r["lineage"], "identity": r["identity"], "aligned_bases": r["aligned_bases"]})
    name = None
    for line in open(fasta):
        if line.startswith(">"):
            name = wanted.get(line[1:].split()[0])
            if name:
                sequences[name] = []
        elif name:
            sequences[name].append(line.strip())

for path in (args.table, args.fasta):
    path.parent.mkdir(parents=True, exist_ok=True)
temporary = args.table.with_name(args.table.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["sample_id", "source", "sequence_id", "lineage", "identity",
                                                "aligned_bases"], delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
os.replace(temporary, args.table)
temporary = args.fasta.with_name(args.fasta.name + ".tmp")
with temporary.open("w") as out:
    for name, parts in sequences.items():
        out.write(f">{name}\n{''.join(parts)}\n")
os.replace(temporary, args.fasta)
print(f"non-host assembled sequences: {len(rows)}")
