#!/usr/bin/env python3
"""Pre-specified contamination test: flowcell concentration of validated presences.

For each non-decoy genome with enough validated presences, the statistic
sum_f n_f^2 / N (n_f = presences on flowcell f) is compared with permutations of
flowcell labels within host species x ocean region strata. Strata with fewer than
the minimum number of flowcells keep their labels fixed; a genome whose presences
lie only in such strata is untestable. BH FDR is applied across tested genomes.
"""
import argparse
import csv
import json
import os
import random
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--analysis", required=True)
parser.add_argument("--grades", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()

rule = json.loads(Path(args.analysis).read_text())["contamination_test"]
libraries = {}
with open(args.manifest, newline="") as handle:
    for r in csv.DictReader(handle):
        if r["include_primary"].lower() != "true":
            continue
        sample = f"{r['study']}__{r['library_id'].split(':')[-1]}"
        flowcells = sorted({"/".join(b.split("/")[:3]) for b in r["sequencing_batches"].split(";") if b})
        libraries[sample] = {"flowcell": ";".join(flowcells), "stratum": (r["species_current"],
                                                                         r["ocean_region"] or "unknown")}
samples = sorted(libraries)
strata = defaultdict(list)
for sample in samples:
    strata[libraries[sample]["stratum"]].append(sample)
permutable = {k for k, members in strata.items()
              if len({libraries[s]["flowcell"] for s in members}) >= rule["min_flowcells_in_stratum"]}

present = defaultdict(set)
taxonomy = {}
with open(args.grades, newline="") as handle:
    for r in csv.DictReader(handle, delimiter="\t"):
        if r["role"] != "decoy" and r["grade"] in ("validated", "high_confidence"):
            present[r["target_id"]].add(r["sample_id"])
            taxonomy[r["target_id"]] = r["taxonomy"]


def concentration(labels, members):
    counts = Counter(labels[s] for s in members)
    return sum(n * n for n in counts.values()) / len(members)


generator = random.Random(rule["seed"])
base = {s: libraries[s]["flowcell"] for s in samples}
permutations = []
for _ in range(rule["permutations"]):
    labels = dict(base)
    for key in permutable:
        members = strata[key]
        shuffled = [base[s] for s in members]
        generator.shuffle(shuffled)
        labels.update(zip(members, shuffled))
    permutations.append(labels)

results = []
for target, members in sorted(present.items()):
    if len(members) < rule["min_validated_presences"]:
        continue
    testable = any(libraries[s]["stratum"] in permutable for s in members)
    observed = concentration(base, members)
    if testable:
        exceed = sum(concentration(labels, members) >= observed for labels in permutations)
        p = (exceed + 1) / (len(permutations) + 1)
    else:
        p = None
    genus = taxonomy[target].split(";")[5][3:] if taxonomy[target] else ""
    results.append({"target_id": target, "taxonomy": taxonomy[target], "validated_presences": len(members),
                    "flowcells": len({base[s] for s in members}), "observed_concentration": round(observed, 4),
                    "p_value": "" if p is None else round(p, 6), "status": "tested" if testable else "untestable",
                    "identity_annotation": str(genus.split("_")[0] in rule["identity_annotation_genera"]).lower()})

tested = sorted((r for r in results if r["status"] == "tested"), key=lambda r: r["p_value"])
m = len(tested)
running = 1.0
for rank, r in reversed(list(enumerate(tested, start=1))):
    running = min(running, r["p_value"] * m / rank)
    r["q_value"] = round(running, 6)
for r in results:
    r.setdefault("q_value", "")
    r["probable_contaminant"] = str(r["status"] == "tested" and r["q_value"] < rule["fdr"]).lower()

fields = ["target_id", "taxonomy", "validated_presences", "flowcells", "observed_concentration", "p_value",
          "q_value", "status", "probable_contaminant", "identity_annotation"]
args.output.parent.mkdir(parents=True, exist_ok=True)
temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(results)
os.replace(temporary, args.output)
print(f"genomes considered={len(results)} tested={m} flagged={sum(r['probable_contaminant'] == 'true' for r in results)}")
