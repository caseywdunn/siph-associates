#!/usr/bin/env python3
"""Check the Phase-3 catalogs against the locked rules and exact accounting."""
import argparse
import csv
import json
import os
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--candidates", required=True)
parser.add_argument("--mag-catalog", required=True)
parser.add_argument("--mag-representatives", required=True, type=Path)
parser.add_argument("--included", required=True)
parser.add_argument("--classification", required=True)
parser.add_argument("--votus", required=True, nargs=2, help="associate and endogenous-candidate vOTU tables")
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


rules = json.loads(Path(args.config).read_text())
errors = []

# MAGs: every retained genome meets the quality rule, is placed in a prokaryotic
# domain, and belongs to exactly one species cluster whose representative exists.
mag_rules = rules["mags"]
candidates = rows(args.candidates)
catalog = rows(args.mag_catalog)
if {m["mag_id"] for m in catalog} != {m["mag_id"] for m in candidates}:
    errors.append("MAG catalog does not account for every quality-selected candidate")
retained = [m for m in catalog if m["catalog_status"] == "retained"]
for mag in catalog:
    if float(mag["completeness"]) < mag_rules["min_completeness"] or \
            float(mag["contamination"]) >= mag_rules["max_contamination_exclusive"]:
        errors.append(f"MAG violates the quality rule: {mag['mag_id']}")
    if not mag["mag_id"].startswith(mag["sample_id"] + "__"):
        errors.append(f"MAG ID does not carry its source library: {mag['mag_id']}")
for mag in retained:
    if mag["gtdb_classification"].split(";")[0] not in {f"d__{d}" for d in mag_rules["require_gtdbtk_domains"]}:
        errors.append(f"retained MAG lacks prokaryotic placement: {mag['mag_id']}")
    if not mag["species_cluster"] or not mag["representative"]:
        errors.append(f"retained MAG has no species cluster: {mag['mag_id']}")
species = {m["species_cluster"]: m["representative"] for m in retained if m["mag_id"] == m["representative"]}
if len(species) != len({m["species_cluster"] for m in retained}):
    errors.append("a species cluster lacks exactly one representative")
for cluster in species:
    if not (args.mag_representatives / f"{cluster}.fa").is_file():
        errors.append(f"missing species representative FASTA: {cluster}")

# Viruses: every included contig meets the inclusion rule, is classified once, and
# appears in exactly one vOTU of its class.
virus_rules = rules["viruses"]
included = rows(args.included)
for virus in included:
    completeness = float(virus["checkv_completeness"]) if virus["checkv_completeness"] not in ("", "NA") else 0.0
    if int(virus["genomad_hallmarks"]) < virus_rules["min_genomad_hallmarks"] or not (
            int(virus["length"]) >= virus_rules["min_length"]
            or completeness >= virus_rules["min_checkv_completeness"]):
        errors.append(f"included virus violates the inclusion rule: {virus['contig_id']}")
classification = {row["contig_id"]: row["virus_class"] for row in rows(args.classification)}
if set(classification) != {v["contig_id"] for v in included}:
    errors.append("viral classification does not match the included set")
seen = Counter()
for path, expected in zip(args.votus, ("associate", "endogenous_candidate")):
    for row in rows(path):
        seen[row["contig_id"]] += 1
        if classification.get(row["contig_id"]) != expected:
            errors.append(f"contig clustered in the wrong class: {row['contig_id']}")
if set(seen) != set(classification) or any(count != 1 for count in seen.values()):
    errors.append("every included viral contig must appear in exactly one vOTU")
votus = {path: {row["votu_id"] for row in rows(path)} for path in args.votus}

lines = [f"status\t{'PASS' if not errors else 'FAIL'}",
         f"mag_candidates\t{len(candidates)}", f"mags_retained\t{len(retained)}",
         f"mag_species\t{len(species)}",
         f"mags_near_complete\t{sum(m['near_complete'] == 'true' for m in retained)}",
         f"viruses_included\t{len(included)}",
         f"viruses_associate\t{sum(c == 'associate' for c in classification.values())}",
         f"viruses_endogenous_candidate\t{sum(c == 'endogenous_candidate' for c in classification.values())}",
         f"votus_associate\t{len(votus[args.votus[0]])}",
         f"votus_endogenous_candidate\t{len(votus[args.votus[1]])}",
         f"errors\t{len(errors)}", *[f"ERROR\t{error}" for error in errors]]
print("\n".join(lines))
if errors:
    raise SystemExit(1)
temporary = args.output.with_name(args.output.name + ".tmp")
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, args.output)
