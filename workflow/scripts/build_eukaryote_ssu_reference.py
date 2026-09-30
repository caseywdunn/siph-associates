#!/usr/bin/env python3
"""Build the competitive SSU reference: SILVA 138.1 NR99 plus this cohort's host SSUs.

SILVA lacks most siphonophore host species, so host reads from conserved 18S
regions otherwise land on arbitrary metazoan references. Every assembled SSU
whose SILVA best hit lies in Siphonophorae is added as a host reference.
"""
import argparse
import csv
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--silva", required=True)
parser.add_argument("--classifications", required=True, nargs="+")
parser.add_argument("--fasta", required=True, type=Path)
parser.add_argument("--taxonomy", required=True, type=Path)
args = parser.parse_args()

hosts = {}
for path in args.classifications:
    assembled = Path(path).with_name(Path(path).name.replace(".ssu_classification.csv", ".assembled_ssu.fasta"))
    with open(path, newline="") as handle:
        classes = {r["OTU"]: r["taxonomy"] for r in csv.DictReader(handle) if "Siphonophorae" in r["taxonomy"]}
    name, sequence = None, []
    for line in open(assembled):
        if line.startswith(">"):
            if name in classes:
                hosts[name] = (classes[name], "".join(sequence))
            name, sequence = line[1:].split()[0].rsplit("_", 1)[0], []
        else:
            sequence.append(line.strip())
    if name in classes:
        hosts[name] = (classes[name], "".join(sequence))

fasta_tmp = args.fasta.with_name(args.fasta.name + ".tmp")
tax_tmp = args.taxonomy.with_name(args.taxonomy.name + ".tmp")
silva = 0
with open(fasta_tmp, "w") as fasta, open(tax_tmp, "w") as taxonomy:
    for line in open(args.silva):
        if line.startswith(">"):
            ref, _, lineage = line[1:].rstrip("\n").partition(" ")
            fasta.write(f">{ref}\n")
            taxonomy.write(f"{ref}\t{lineage}\n")
            silva += 1
        else:
            fasta.write(line)
    for otu, (lineage, sequence) in sorted(hosts.items()):
        ranks = lineage.split(";")
        host = ";".join(ranks[:ranks.index("Siphonophorae") + 1] + ["host siphonophore", otu])
        ref = f"HOST|{otu}"
        fasta.write(f">{ref}\n{sequence}\n")
        taxonomy.write(f"{ref}\t{host}\n")
os.replace(fasta_tmp, args.fasta)
os.replace(tax_tmp, args.taxonomy)
print(f"references: SILVA={silva} host={len(hosts)}")
