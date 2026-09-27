#!/usr/bin/env python3
"""Apply the locked viral inclusion rule to per-library geNomad and CheckV outputs."""
import argparse
import csv
import json
import os
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--membership", required=True)
parser.add_argument("--results", required=True, type=Path)
parser.add_argument("--config", required=True)
parser.add_argument("--fasta", required=True, type=Path)
parser.add_argument("--table", required=True, type=Path)
args = parser.parse_args()

rules = json.loads(Path(args.config).read_text())["viruses"]
fields = ["contig_id", "sample_id", "length", "genomad_hallmarks", "genomad_score", "genomad_topology",
          "checkv_quality", "checkv_completeness", "checkv_completeness_method", "taxonomy", "phylum"]


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


with open(args.membership, newline="") as handle:
    samples = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")
               if row["assembly_eligible"] == "true"]
included = []
fasta_tmp = args.fasta.with_name(args.fasta.name + ".tmp")
with fasta_tmp.open("w") as out:
    for sample in samples:
        checkv = {row["contig_id"]: row for row in rows(args.results / "checkv" / f"{sample}.tsv")}
        keep = set()
        for row in rows(args.results / "viruses" / f"{sample}.tsv"):
            if not row.get("length"):
                continue
            quality = checkv[row["seq_name"]]
            completeness = float(quality["completeness"]) if quality["completeness"] not in ("", "NA") else 0.0
            long_enough = int(row["length"]) >= rules["min_length"]
            complete_enough = completeness >= rules["min_checkv_completeness"]
            # Terminal repeats never qualify alone: a hallmark gene is always required.
            if int(row["n_hallmarks"]) < rules["min_genomad_hallmarks"] or not (long_enough or complete_enough):
                continue
            ranks = row["taxonomy"].split(";")
            included.append({
                "contig_id": row["seq_name"], "sample_id": sample, "length": row["length"],
                "genomad_hallmarks": row["n_hallmarks"], "genomad_score": row["virus_score"],
                "genomad_topology": row["topology"], "checkv_quality": quality["checkv_quality"],
                "checkv_completeness": quality["completeness"],
                "checkv_completeness_method": quality["completeness_method"],
                "taxonomy": row["taxonomy"], "phylum": ranks[3] if len(ranks) > 3 and ranks[3] else "Unclassified",
            })
            keep.add(row["seq_name"])
        write = False
        for line in open(args.results / "viruses" / f"{sample}.fna"):
            if line.startswith(">"):
                write = line[1:].split()[0] in keep
            if write:
                out.write(line)
os.replace(fasta_tmp, args.fasta)

temporary = args.table.with_name(args.table.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(included)
os.replace(temporary, args.table)
print(f"included {len(included)} viral contigs from {len({row['sample_id'] for row in included})} libraries")
