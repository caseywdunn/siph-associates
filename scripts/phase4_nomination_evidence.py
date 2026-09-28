#!/usr/bin/env python3
"""Tabulate the evidence behind the Phase-4 bacterial nomination decision.

Reads the accepted Phase-2 cohort nominations and the Phase-3 MAG catalog and
writes the tables cited in docs/phase4_catalog_decisions.md to
data/results/phase4_catalog/nomination_evidence/.
"""
from __future__ import annotations

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "data" / "results"
OUT = RESULTS / "phase4_catalog" / "nomination_evidence"
TAXONOMY = Path("/gpfs/ycga/work/dunn/cwd7/databases/gtdbtk_r220/release220/taxonomy/gtdb_taxonomy.tsv")
CONTAMINANT_GENERA = {"Cutibacterium", "Pelomonas", "Bradyrhizobium", "Mesorhizobium", "Ralstonia",
                      "Staphylococcus", "Corynebacterium", "Streptococcus", "Escherichia", "Acinetobacter",
                      "Delftia", "Sphingomonas", "Methylobacterium", "Micrococcus"}


def genus_name(lineage):
    """GTDB genus without its alphabetic suffix (Vibrio_A -> Vibrio) for cross-taxonomy matching."""
    genus = lineage.split(";")[5][3:] if lineage else ""
    return re.sub(r"_[A-Z]+$", "", genus)


def write(name, header, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n# {name}")
    for row in [header, *rows]:
        print("\t".join(str(value) for value in row))


taxonomy = {}
for line in TAXONOMY.open():
    accession, lineage = line.rstrip("\n").split("\t")
    taxonomy[accession[3:]] = lineage

sylph, bracken, phyloflash = [], defaultdict(dict), defaultdict(Counter)
with open(RESULTS / "aggregation" / "cohort" / "candidate_nominations.tsv", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["source"] == "sylph":
            sylph.append(row)
        elif row["source"] == "Bracken" and row["rank"] == "genus":
            bracken[row["sample_id"]][row["candidate_name"]] = float(row["estimated_reads"] or 0)
        elif row["source"] == "phyloFlash":
            ranks = row["candidate_id"].split(";")
            if ranks[0] in ("Bacteria", "Archaea") and len(ranks) >= 6:
                phyloflash[row["sample_id"]][ranks[5]] += float(row["raw_support"] or 0)

hits = []
for row in sylph:
    accession = re.sub(r"_genomic\.fna\.gz$", "", row["candidate_name"])
    source = "OceanDNA" if accession.startswith("OceanDNA") else "GTDB"
    lineage = taxonomy.get(accession, "")
    genus = genus_name(lineage)
    hits.append({
        "sample": row["sample_id"], "accession": accession, "source": source, "genus": genus,
        "species": lineage.split(";")[6] if lineage else "", "ani": float(row["adjusted_ani"]),
        "bracken": bracken[row["sample_id"]].get(genus, 0) if genus else 0,
        "phyloflash": phyloflash[row["sample_id"]].get(genus, 0) if genus else 0,
    })
corroborated = lambda h: h["bracken"] >= 10 or h["phyloflash"] >= 1
median = lambda values: sorted(values)[len(values) // 2] if values else ""

gtdb = [h for h in hits if h["source"] == "GTDB"]
confirmed = {h["accession"] for h in gtdb if corroborated(h)}
write("sylph_nominations.tsv", ["source", "hits", "genomes", "libraries", "median_ani"],
      [[s, len(group), len({h["accession"] for h in group}), len({h["sample"] for h in group}),
        median([h["ani"] for h in group])]
       for s, group in (("GTDB", gtdb), ("OceanDNA", [h for h in hits if h["source"] == "OceanDNA"]),
                        ("all", hits))])
write("gtdb_corroboration.tsv", ["evidence", "hits_corroborated", "gtdb_hits"],
      [["Bracken genus >=10 reads", sum(h["bracken"] >= 10 for h in gtdb), len(gtdb)],
       ["phyloFlash genus >=1 read", sum(h["phyloflash"] >= 1 for h in gtdb), len(gtdb)],
       ["either", sum(corroborated(h) for h in gtdb), len(gtdb)]])
unconfirmed = [h for h in gtdb if h["accession"] not in confirmed]
write("gtdb_uncorroborated.tsv", ["genomes", "median_ani", "corroborated_genomes", "corroborated_median_ani",
                                  "top_genera"],
      [[len({h["accession"] for h in unconfirmed}), median([h["ani"] for h in unconfirmed]), len(confirmed),
        median([h["ani"] for h in gtdb if corroborated(h)]),
        ",".join(f"{g}:{n}" for g, n in Counter(h["genus"] for h in unconfirmed).most_common(10))]])
contaminants = [h for h in gtdb if h["genus"] in CONTAMINANT_GENERA]
write("contaminant_genera.tsv", ["hits", "genomes", "libraries", "genera"],
      [[len(contaminants), len({h["accession"] for h in contaminants}), len({h["sample"] for h in contaminants}),
        ",".join(f"{g}:{n}" for g, n in Counter(h["genus"] for h in contaminants).most_common())]])

with open(RESULTS / "phase3_catalog" / "mags" / "mag_catalog.tsv", newline="") as handle:
    mags = [row for row in csv.DictReader(handle, delimiter="\t")
            if row["catalog_status"] == "retained" and row["mag_id"] == row["representative"]]
named = [m["gtdb_classification"].split(";")[6] for m in mags if m["gtdb_classification"].split(";")[6] != "s__"]
nominated = {h["species"] for h in gtdb}
write("mag_overlap.tsv", ["mag_species", "named_species", "named_also_nominated_by_sylph", "species"],
      [[len(mags), len(named), sum(s in nominated for s in named), ",".join(named)]])
