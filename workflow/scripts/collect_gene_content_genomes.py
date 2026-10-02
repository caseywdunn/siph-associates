#!/usr/bin/env python3
"""Gather the phylogeny's Mycoplasmatales genomes and every GTDB r220 species representative of one family.

GTDB genomes are named GTDB_<accession>; an accession already present as an external genome is not repeated.
"""
import argparse
import csv
import gzip
import os
import shutil
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--tree-genomes", required=True, type=Path)
parser.add_argument("--labels", required=True)
parser.add_argument("--gtdb", required=True, type=Path, help="GTDB-Tk r220 data directory")
parser.add_argument("--taxonomy", required=True, help="taxonomy table relative to --gtdb")
parser.add_argument("--family", required=True)
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--table", required=True, type=Path)
args = parser.parse_args()

with open(args.labels, newline="") as handle:
    labels = [r for r in csv.DictReader(handle, delimiter="\t") if r["status"] == "included"]
with open(args.gtdb / "skani" / "genome_paths.tsv") as handle:
    paths = {line.split()[0]: args.gtdb / "skani" / line.split()[1] / line.split()[0] for line in handle if line.strip()}

args.genomes.mkdir(parents=True)
records = []
for label in labels:
    shutil.copyfile(args.tree_genomes / f"{label['genome']}.fa", args.genomes / f"{label['genome']}.fa")
    records.append({"genome": label["genome"], "group": label["group"], "host_or_source": label["host_or_source"],
                    "classification": label["classification"]})
external = {r["genome"].removeprefix("EXT_") for r in records}

for line in open(args.gtdb / args.taxonomy):
    accession, taxonomy = line.rstrip("\n").split("\t")
    if f"f__{args.family}" not in taxonomy.split(";"):
        continue
    accession = accession.removeprefix("RS_").removeprefix("GB_")
    if accession in external:
        continue
    with gzip.open(paths[f"{accession}_genomic.fna.gz"], "rt") as source, \
            open(args.genomes / f"GTDB_{accession}.fa", "w") as out:
        shutil.copyfileobj(source, out)
    records.append({"genome": f"GTDB_{accession}", "group": "GTDB r220 species representative",
                    "host_or_source": taxonomy.split(";")[-1].removeprefix("s__"), "classification": taxonomy})

temporary = args.table.with_name(args.table.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
os.replace(temporary, args.table)
print(f"genomes={len(records)} (tree={len(labels)}, gtdb={len(records) - len(labels)})")
