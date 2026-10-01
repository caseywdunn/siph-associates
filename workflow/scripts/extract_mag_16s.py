#!/usr/bin/env python3
"""Write the 16S rRNA genes barrnap found on contigs of the selected MAGs as FASTA.

Headers are <species_cluster>|<mag_id>|<contig>:<start>-<end><strand>.
"""
import argparse
import csv
import os
import tarfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--mag-catalog", required=True)
parser.add_argument("--results", required=True, type=Path, help="phase3_cohort directory")
parser.add_argument("--order", required=True)
parser.add_argument("--min-length", type=int, default=1200)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()

complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
with open(args.mag_catalog, newline="") as handle:
    mags = [r for r in csv.DictReader(handle, delimiter="\t")
            if r["catalog_status"] == "retained" and f"o__{args.order}" in r["gtdb_classification"]]
written = 0
temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w") as out:
    for mag in mags:
        with tarfile.open(args.results / "bins" / f"{mag['sample_id']}.tar.gz") as archive:
            text = archive.extractfile(f"{mag['mag_id']}.fa").read().decode()
        contigs, name = {}, None
        for line in text.splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                contigs[name] = []
            else:
                contigs[name].append(line.strip())
        for line in open(args.results / "markers" / f"{mag['sample_id']}.gff"):
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in contigs or "Name=16S_rRNA" not in f[8] or "kingdom=bac" not in f[8]:
                continue
            start, end = int(f[3]), int(f[4])
            if end - start + 1 < args.min_length:
                continue
            gene = "".join(contigs[f[0]])[start - 1:end]
            if f[6] == "-":
                gene = gene.translate(complement)[::-1]
            out.write(f">{mag['species_cluster']}|{mag['mag_id']}|{f[0]}:{start}-{end}{f[6]}\n{gene}\n")
            written += 1
os.replace(temporary, args.output)
print(f"16S genes >= {args.min_length} bp: {written}")
