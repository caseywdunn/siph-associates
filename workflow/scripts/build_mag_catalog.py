#!/usr/bin/env python3
"""Keep prokaryotic MAGs, cluster species by ANI, and write representatives."""
import argparse
import csv
import json
import os
import shutil
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--candidates", required=True)
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--gtdbtk", required=True, nargs="+")
parser.add_argument("--ani", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--catalog", required=True, type=Path)
parser.add_argument("--representatives", required=True, type=Path)
args = parser.parse_args()

rules = json.loads(Path(args.config).read_text())["mags"]
domains = {f"d__{domain}" for domain in rules["require_gtdbtk_domains"]}
ani_min, af_min = rules["dereplication"]["ani_min"], rules["dereplication"]["align_fraction_min"]


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


mags = rows(args.candidates)
classification = {}
for path in args.gtdbtk:
    for row in rows(path):
        classification[row["user_genome"]] = row
for mag in mags:
    record = classification.get(mag["mag_id"], {})
    mag["gtdb_classification"] = record.get("classification", "not_placed")
    mag["gtdb_closest_ani"] = record.get("closest_genome_ani", "")
    mag["gtdb_red"] = record.get("red_value", "")
    placed = mag["gtdb_classification"].split(";")[0] in domains
    mag["catalog_status"] = "retained" if placed else "excluded_no_prokaryotic_placement"

# Pairwise ANI among candidates; species links also need the configured alignment-fraction rule.
linked = defaultdict(set)
for row in rows(args.ani):
    a, b = Path(row["Ref_file"]).name[:-3], Path(row["Query_file"]).name[:-3]
    fractions = (float(row["Align_fraction_ref"]), float(row["Align_fraction_query"]))
    # "at_least_one": an incomplete MAG cannot cover its partner, so require coverage of one genome only.
    covered = max(fractions) if rules["dereplication"].get("align_fraction_rule") == "at_least_one" else min(fractions)
    if float(row["ANI"]) >= ani_min and covered >= af_min:
        linked[a].add(b)
        linked[b].add(a)

# Greedy centroid clustering in descending quality order: every member meets the
# species threshold with its representative, and the representative has the best score.
retained = sorted((m for m in mags if m["catalog_status"] == "retained"),
                  key=lambda m: (-float(m["quality_score"]), m["mag_id"]))
representatives = []
for mag in retained:
    representative = next((r for r in representatives if r["mag_id"] in linked[mag["mag_id"]]), None)
    if representative is None:
        representatives.append(mag)
        representative = mag
    mag["species_cluster"] = f"MAGSP{representatives.index(representative) + 1:04d}"
    mag["representative"] = representative["mag_id"]
for mag in mags:
    mag.setdefault("species_cluster", "")
    mag.setdefault("representative", "")

args.representatives.mkdir(parents=True)
for representative in representatives:
    shutil.copyfile(args.genomes / f"{representative['mag_id']}.fa",
                    args.representatives / f"{representative['species_cluster']}.fa")

fields = list(mags[0]) if mags else ["mag_id"]
temporary = args.catalog.with_name(args.catalog.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(sorted(mags, key=lambda m: (m["species_cluster"] or "~", m["mag_id"])))
os.replace(temporary, args.catalog)

# Single-linkage count for comparison with the decision-record estimate.
parent = {m["mag_id"]: m["mag_id"] for m in retained}
def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x
for a, targets in linked.items():
    for b in targets:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
print(f"candidates={len(mags)} retained={len(retained)} species={len(representatives)} "
      f"single_linkage_species={len({find(x) for x in parent})}")
