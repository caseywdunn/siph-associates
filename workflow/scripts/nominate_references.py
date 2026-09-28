#!/usr/bin/env python3
"""Nominate reference genomes from cohort sylph hits and draw decoy candidates."""
import argparse
import csv
import json
import os
import random
import re
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--nominations", required=True)
parser.add_argument("--taxonomy", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--references", required=True, type=Path)
parser.add_argument("--decoys", required=True, type=Path)
args = parser.parse_args()

settings = json.loads(Path(args.config).read_text())
rules, decoy_rules = settings["nomination"], settings["decoys"]
taxonomy = {}
for line in open(args.taxonomy):
    accession, lineage = line.rstrip("\n").split("\t")
    taxonomy[accession[3:]] = lineage


def genus_name(lineage):
    genus = lineage.split(";")[5][3:] if lineage else ""
    return re.sub(r"_[A-Z]+$", "", genus)


hits, bracken, phyloflash = [], defaultdict(dict), defaultdict(Counter)
with open(args.nominations, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["source"] == "sylph":
            hits.append(row)
        elif row["source"] == "Bracken" and row["rank"] == "genus":
            bracken[row["sample_id"]][row["candidate_name"]] = float(row["estimated_reads"] or 0)
        elif row["source"] == "phyloFlash":
            ranks = row["candidate_id"].split(";")
            if ranks[0] in ("Bacteria", "Archaea") and len(ranks) >= 6:
                phyloflash[row["sample_id"]][ranks[5]] += float(row["raw_support"] or 0)

genomes = defaultdict(list)
for row in hits:
    if float(row["adjusted_ani"]) < rules["sylph_min_adjusted_ani"]:
        continue
    name = row["candidate_name"]
    genome_id = name.removesuffix(".fa.gz") if name.startswith("OceanDNA") else re.sub(r"_genomic\.fna\.gz$", "", name)
    genomes[genome_id].append(row)

fields = ["genome_id", "source", "gtdb_taxonomy", "libraries", "max_ani", "median_ani",
          "bracken_corroborated_libraries", "phyloflash_corroborated_libraries", "library_list"]
rows = []
for genome_id, group in sorted(genomes.items()):
    lineage = taxonomy.get(genome_id, "")
    genus = genus_name(lineage)
    anis = sorted(float(r["adjusted_ani"]) for r in group)
    samples = sorted({r["sample_id"] for r in group})
    rows.append({
        "genome_id": genome_id,
        "source": "OceanDNA" if genome_id.startswith("OceanDNA") else "GTDB",
        "gtdb_taxonomy": lineage, "libraries": len(samples), "max_ani": anis[-1],
        "median_ani": anis[len(anis) // 2],
        "bracken_corroborated_libraries": sum(
            bracken[s].get(genus, 0) >= rules["evidence_bracken_genus_min_reads"] for s in samples) if genus else 0,
        "phyloflash_corroborated_libraries": sum(
            phyloflash[s].get(genus, 0) >= rules["evidence_phyloflash_genus_min_reads"] for s in samples) if genus else 0,
        "library_list": ",".join(samples),
    })


def write(path, header, records):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    os.replace(temporary, path)


write(args.references, fields, rows)

# Decoy candidates: a seeded ordering of GTDB representatives within each family,
# oversampled threefold so later QC and ANI exclusions can be replaced in order.
generator = random.Random(decoy_rules["seed"])
nominated = set(genomes)
decoys = []
for family, quota in decoy_rules["families"].items():
    pool = sorted(a for a, lineage in taxonomy.items() if lineage.split(";")[4] == family and a not in nominated)
    for rank, accession in enumerate(generator.sample(pool, min(len(pool), 3 * quota)), start=1):
        decoys.append({"genome_id": accession, "family": family, "draw_order": rank, "quota": quota,
                       "gtdb_taxonomy": taxonomy[accession]})
write(args.decoys, ["genome_id", "family", "draw_order", "quota", "gtdb_taxonomy"], decoys)
print(f"nominated {len(rows)} genomes ({Counter(r['source'] for r in rows)}); decoy candidates {len(decoys)}")
