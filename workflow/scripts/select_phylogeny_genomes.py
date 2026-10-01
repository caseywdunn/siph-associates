#!/usr/bin/env python3
"""Choose genomes for the Mycoplasmatales tree and write the GTDB-Tk input directory and label table.

Cohort genomes: every retained Phase-3 MAG species placed in Mycoplasmatales.
External genomes: those GTDB-Tk places in Mycoplasmatales with CheckM2
completeness at or above the minimum. Every external genome is accounted for
in the label table, with its inclusion status.
"""
import argparse
import csv
import os
import shutil
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--external", required=True, help="config/mycoplasmatales_external.tsv")
parser.add_argument("--external-genomes", required=True, type=Path)
parser.add_argument("--checkm2", required=True)
parser.add_argument("--gtdbtk", required=True, nargs="+")
parser.add_argument("--mag-catalog", required=True)
parser.add_argument("--mag-representatives", required=True, type=Path)
parser.add_argument("--order", required=True)
parser.add_argument("--min-completeness", required=True, type=float)
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--labels", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


quality = {r["Name"]: r for r in rows(args.checkm2)}
placement = {r["user_genome"]: r["classification"] for path in args.gtdbtk for r in rows(path)}
args.genomes.mkdir(parents=True)
labels = []
for mag in rows(args.mag_catalog):
    if (mag["catalog_status"] == "retained" and mag["mag_id"] == mag["representative"]
            and f"o__{args.order}" in mag["gtdb_classification"]):
        shutil.copyfile(args.mag_representatives / f"{mag['species_cluster']}.fa",
                        args.genomes / f"{mag['species_cluster']}.fa")
        labels.append({"genome": mag["species_cluster"], "group": "siphonophore (this study)",
                       "host_or_source": f"source library {mag['sample_id']}", "classification": mag["gtdb_classification"],
                       "completeness": mag["completeness"], "contamination": mag["contamination"], "status": "included"})
for ext in rows(args.external):
    accession = ext["accession"]
    name = f"EXT_{accession}"
    q = quality.get(name, {})
    completeness = float(q.get("Completeness", 0) or 0)
    classification = placement.get(name, "not_placed")
    if f"o__{args.order}" not in classification:
        status = "excluded_outside_order"
    elif completeness < args.min_completeness:
        status = "excluded_low_completeness"
    else:
        status = "included"
        shutil.copyfile(args.external_genomes / f"{name}.fa", args.genomes / f"{name}.fa")
    labels.append({"genome": name, "group": ext["reason"], "host_or_source": ext["host_or_source"],
                   "classification": classification, "completeness": q.get("Completeness", ""),
                   "contamination": q.get("Contamination", ""), "status": status})

temporary = args.labels.with_name(args.labels.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(labels[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(labels)
os.replace(temporary, args.labels)
print(f"included {sum(l['status'] == 'included' for l in labels)} of {len(labels)} genomes")
