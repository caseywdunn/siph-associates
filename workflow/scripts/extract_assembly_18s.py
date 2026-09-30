#!/usr/bin/env python3
"""Write the 18S rRNA genes barrnap found on a Phase-3 assembly as FASTA.

Only eukaryotic 18S features of at least --min-length bp are kept; shorter
partial genes come from conserved regions that cannot identify a lineage.
"""
import argparse
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--gff", required=True)
parser.add_argument("--contigs", required=True)
parser.add_argument("--min-length", type=int, required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()

features = []
for line in open(args.gff):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 9 and "Name=18S_rRNA" in f[8] and "kingdom=euk" in f[8]:
        start, end = int(f[3]), int(f[4])
        if end - start + 1 >= args.min_length:
            features.append((f[0], start, end, f[6]))
wanted = {contig for contig, *_ in features}
sequences, name, chunks = {}, None, []
for line in open(args.contigs):
    if line.startswith(">"):
        if name in wanted:
            sequences[name] = "".join(chunks)
        name, chunks = line[1:].split()[0], []
    elif name in wanted:
        chunks.append(line.strip())
if name in wanted:
    sequences[name] = "".join(chunks)

complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w") as out:
    for i, (contig, start, end, strand) in enumerate(features, start=1):
        gene = sequences[contig][start - 1:end]
        if strand == "-":
            gene = gene.translate(complement)[::-1]
        out.write(f">{contig}:{start}-{end}{strand}\n{gene}\n")
os.replace(temporary, args.output)
print(f"18S genes >= {args.min_length} bp: {len(features)}")
