#!/usr/bin/env python3
"""Check the frozen Phase-4 catalogs: attribution, rules, checksums, and exact accounting."""
import argparse
import csv
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--references", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--mag-catalog", required=True)
parser.add_argument("--viral-manifests", required=True, nargs=2)
parser.add_argument("--viral-fastas", required=True, nargs=2)
parser.add_argument("--links", required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


settings = json.loads(Path(args.config).read_text())
qc = settings["reference_qc"]
errors = []
manifest = rows(args.manifest)
members = [m for m in manifest if m["status"] == "representative"]
real = [m for m in members if m["role"] != "decoy"]
decoys = [m for m in members if m["role"] == "decoy"]

# Every nominated genome is accounted for exactly once.
nominated = {r["genome_id"] for r in rows(args.references)}
accounted = Counter(m["genome_id"] for m in manifest if m["role"] == "reference")
if set(accounted) != nominated or any(c != 1 for c in accounted.values()):
    errors.append("manifest does not account for every nominated reference exactly once")
mag_species = {r["species_cluster"] for r in rows(args.mag_catalog)
               if r["catalog_status"] == "retained" and r["mag_id"] == r["representative"]}
if {m["genome_id"] for m in manifest if m["role"] == "mag"} != mag_species:
    errors.append("manifest does not account for every MAG species")

# Members meet the rules, have unique IDs, and match the FASTA contig prefixes.
ids = [m["catalog_id"] for m in members]
if len(ids) != len(set(ids)) or any(not i for i in ids):
    errors.append("catalog IDs are missing or duplicated")
for member in members:
    if float(member["completeness"]) < qc["min_completeness"] or \
            float(member["contamination"]) >= qc["max_contamination_exclusive"]:
        errors.append(f"catalog genome violates the quality rule: {member['catalog_id']}")
clustered = [m for m in manifest if m["status"].startswith("clustered_into:")]
listed = Counter(g for m in real for g in m["cluster_members"].split(",") if g)
for entry in clustered:
    if listed[entry["genome_id"]] != 1:
        errors.append(f"clustered genome is not listed under exactly one representative: {entry['genome_id']}")
quota = settings["decoys"]["families"]
families = Counter(m["evidence"].split(";")[0].split("=")[1] for m in decoys)
if families != Counter(quota):
    errors.append(f"decoy family counts {dict(families)} differ from the locked quotas")
contigs = Counter()
for line in open(args.fasta):
    if line.startswith(">"):
        contigs[line[1:].split("|")[0]] += 1
for member in members:
    if contigs[member["catalog_id"]] != int(member["contigs"]):
        errors.append(f"FASTA contig count differs from the manifest: {member['catalog_id']}")
if set(contigs) != set(ids):
    errors.append("catalog FASTA contains genomes absent from the manifest")

# Viral catalogs: checksums and disjoint IDs.
viral_ids = []
for manifest_path, fasta_path in zip(args.viral_manifests, args.viral_fastas):
    entries = {r["votu_id"]: r for r in rows(manifest_path)}
    viral_ids.extend(entries)
    sequences, name = {}, None
    for line in open(fasta_path):
        if line.startswith(">"):
            name = line[1:].strip()
            sequences[name] = []
        else:
            sequences[name].append(line.strip())
    for votu, parts in sequences.items():
        if hashlib.sha256("".join(parts).encode()).hexdigest() != entries.get(votu, {}).get("sequence_sha256"):
            errors.append(f"vOTU sequence checksum mismatch: {votu}")
    if set(sequences) != set(entries):
        errors.append(f"viral FASTA and manifest disagree: {manifest_path}")
if len(viral_ids) != len(set(viral_ids)):
    errors.append("viral catalog IDs overlap between classes")
links = rows(args.links)
if any(l["votu_id"] not in set(viral_ids) or l["host_catalog_id"] not in contigs for l in links):
    errors.append("CRISPR links reference unknown vOTUs or hosts")

lines = [f"status\t{'PASS' if not errors else 'FAIL'}",
         f"catalog_version\t{settings['catalog_version']}",
         f"nominated_references\t{len(nominated)}",
         f"catalog_genomes\t{len(real)}",
         f"catalog_mags\t{sum(m['role'] == 'mag' for m in real)}",
         f"catalog_references\t{sum(m['role'] == 'reference' for m in real)}",
         f"decoys\t{len(decoys)}",
         f"references_excluded\t{sum(m['role'] == 'reference' and not m['status'].startswith(('representative', 'clustered')) for m in manifest)}",
         f"references_clustered\t{len(clustered)}",
         f"associate_votus\t{len(rows(args.viral_manifests[0]))}",
         f"endogenous_candidate_votus\t{len(rows(args.viral_manifests[1]))}",
         f"crispr_host_links\t{len(links)}",
         f"errors\t{len(errors)}", *[f"ERROR\t{e}" for e in errors]]
print("\n".join(lines))
if errors:
    raise SystemExit(1)
temporary = args.output.with_name(args.output.name + ".tmp")
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, args.output)
