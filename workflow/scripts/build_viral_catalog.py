#!/usr/bin/env python3
"""Freeze a Phase-3 vOTU representative set as a versioned viral catalog with checksums."""
import argparse
import csv
import hashlib
import os
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--representatives", required=True)
parser.add_argument("--votus", required=True)
parser.add_argument("--catalog-class", required=True)
parser.add_argument("--fasta", required=True, type=Path)
parser.add_argument("--manifest", required=True, type=Path)
args = parser.parse_args()

with open(args.votus, newline="") as handle:
    members = defaultdict(list)
    for row in csv.DictReader(handle, delimiter="\t"):
        members[row["votu_id"]].append(row)

sequences, name = {}, None
for line in open(args.representatives):
    if line.startswith(">"):
        name = line[1:].split()[0]
        sequences[name] = []
    else:
        sequences[name].append(line.strip())
if set(sequences) != set(members):
    raise SystemExit("representative FASTA and vOTU table disagree")

fasta_tmp = args.fasta.with_name(args.fasta.name + ".tmp")
manifest_tmp = args.manifest.with_name(args.manifest.name + ".tmp")
with fasta_tmp.open("w") as fasta, manifest_tmp.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["votu_id", "catalog_class", "representative_contig", "length", "members", "libraries",
                     "phylum", "taxonomy", "checkv_quality", "endogenous_reasons", "sequence_sha256"])
    for votu in sorted(sequences):
        sequence = "".join(sequences[votu])
        rep = next(r for r in members[votu] if r["is_representative"] == "true")
        fasta.write(f">{votu}\n")
        for i in range(0, len(sequence), 80):
            fasta.write(sequence[i:i + 80] + "\n")
        writer.writerow([votu, args.catalog_class, rep["contig_id"], len(sequence), len(members[votu]),
                         len({r["sample_id"] for r in members[votu]}), rep["phylum"], rep["taxonomy"],
                         rep["checkv_quality"], rep["endogenous_reasons"],
                         hashlib.sha256(sequence.encode()).hexdigest()])
os.replace(fasta_tmp, args.fasta)
os.replace(manifest_tmp, args.manifest)
print(f"{args.catalog_class}: {len(sequences)} vOTUs")
