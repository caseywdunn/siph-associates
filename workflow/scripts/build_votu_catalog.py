#!/usr/bin/env python3
"""Name CheckV aniclust clusters as vOTUs and write their representative sequences."""
import argparse
import csv
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--clusters", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--classification", required=True)
parser.add_argument("--prefix", required=True)
parser.add_argument("--catalog", required=True, type=Path)
parser.add_argument("--representatives", required=True, type=Path)
args = parser.parse_args()

with open(args.classification, newline="") as handle:
    viruses = {row["contig_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
clusters = []
for line in open(args.clusters):
    if line.strip():
        representative, members = line.rstrip("\n").split("\t")
        clusters.append((representative, members.split(",")))
# Stable IDs: order vOTUs by representative length, then name.
clusters.sort(key=lambda c: (-int(viruses[c[0]]["length"]), c[0]))
names = {representative: f"{args.prefix}{i:05d}" for i, (representative, _) in enumerate(clusters, start=1)}

rows = []
for representative, members in clusters:
    for member in members:
        virus = viruses[member]
        rows.append({"votu_id": names[representative], "contig_id": member,
                     "is_representative": str(member == representative).lower(),
                     "cluster_size": len(members),
                     "cluster_libraries": len({viruses[m]["sample_id"] for m in members}),
                     **{key: virus[key] for key in ("sample_id", "length", "checkv_quality", "checkv_completeness",
                                                    "genomad_hallmarks", "phylum", "taxonomy",
                                                    "host_reference_tested", "host_aligned_fraction",
                                                    "endogenous_reasons")}})
fields = list(rows[0]) if rows else ["votu_id", "contig_id"]
temporary = args.catalog.with_name(args.catalog.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
os.replace(temporary, args.catalog)

temporary = args.representatives.with_name(args.representatives.name + ".tmp")
with temporary.open("w") as out:
    write = False
    for line in open(args.fasta):
        if line.startswith(">"):
            contig = line[1:].split()[0]
            write = contig in names
            if write:
                line = f">{names[contig]} {contig}\n"
        if write:
            out.write(line)
os.replace(temporary, args.representatives)
print(f"{args.prefix}: {len(clusters)} clusters from {len(rows)} contigs")
