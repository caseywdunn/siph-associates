#!/usr/bin/env python3
"""Assign pre-specified evidence grades to every library x genome and library x vOTU pair."""
import argparse
import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--presence", required=True)
parser.add_argument("--analysis", required=True)
parser.add_argument("--bacterial", required=True)
parser.add_argument("--viral", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--nominated", required=True)
parser.add_argument("--mags", required=True)
parser.add_argument("--votus", required=True, nargs="+")
parser.add_argument("--support-dir", required=True, type=Path)
parser.add_argument("--out-bacterial", required=True, type=Path)
parser.add_argument("--out-viral", required=True, type=Path)
args = parser.parse_args()

presence = json.loads(Path(args.presence).read_text())
support_rule = json.loads(Path(args.analysis).read_text())["evidence_grades"]["high_confidence_assembly_support"]
bac, vir = presence["bacteria"], presence["viruses"]
READ_LENGTH = 150


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, fields, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    os.replace(temporary, path)


# Catalog genome for every member genome (references and MAG species).
manifest = rows(args.manifest)
catalog_of = {}
for entry in manifest:
    if entry["status"] == "representative":
        for member in entry["cluster_members"].split(","):
            catalog_of[member] = entry["catalog_id"]

# Phase-2 nomination: libraries where sylph nominated a member of the catalog genome.
nominated = defaultdict(set)
for row in rows(args.nominated):
    target = catalog_of.get(row["genome_id"])
    if target:
        for sample in row["library_list"].split(","):
            nominated[target].add(sample)

# MAG species assembled from each library.
mag_libraries = defaultdict(set)
for mag in rows(args.mags):
    if mag["catalog_status"] == "retained":
        target = catalog_of.get(mag["species_cluster"])
        if target:
            mag_libraries[target].add(mag["sample_id"])


def union_length(intervals):
    total, end = 0, -1
    for start, stop in sorted(intervals):
        if stop <= end:
            continue
        total += stop - max(start, end)
        end = stop
    return total


# Same-library assembly support: union of catalog bases aligned at >= the locked identity.
support = {}
for paf in args.support_dir.glob("*.paf"):
    sample = paf.stem
    spans = defaultdict(list)
    for line in open(paf):
        f = line.split("\t")
        matches, length = int(f[9]), int(f[10])
        if length and 100 * matches / length >= support_rule["min_identity_percent"]:
            spans[(f[5].split("|")[0], f[5])].append((int(f[7]), int(f[8])))
    per_genome = defaultdict(int)
    for (genome, _), intervals in spans.items():
        per_genome[genome] += union_length(intervals)
    support[sample] = per_genome


def bacterial_grade(breadth, reads, supported, threshold):
    if breadth >= threshold and reads >= bac["validated"]["min_reads"]:
        return "high_confidence" if supported else "validated"
    return "trace" if breadth >= bac["trace"]["min_breadth"] else "none"


records = []
for row in rows(args.bacterial):
    sample, target = row["sample_id"], row["target_id"]
    breadth, reads, length = float(row["covered_fraction"]), int(row["read_count"]), int(row["length"])
    depth = reads * READ_LENGTH / length if length else 0.0
    expected = 1 - math.exp(-depth)
    ratio = breadth / expected if expected > 0 else 0.0
    aligned = support.get(sample, {}).get(target, 0)
    supported = sample in mag_libraries[target] or aligned >= support_rule["min_aligned_bases"]
    grade = bacterial_grade(breadth, reads, supported, bac["validated"]["min_breadth"])
    records.append({
        "sample_id": sample, "target_id": target, "role": row["role"], "taxonomy": row["taxonomy"],
        "breadth": round(breadth, 6), "reads": reads, "read_depth": round(depth, 4),
        "breadth_ratio": round(ratio, 4),
        "divergent_strain": str(grade in ("validated", "high_confidence")
                                and ratio < bac["divergent_strain_annotation"]["flag_below"]).lower(),
        "nominated": str(sample in nominated[target]).lower(),
        "mag_from_library": str(sample in mag_libraries[target]).lower(),
        "assembly_support_bases": aligned, "grade": grade,
        **{f"grade_at_{int(t * 100)}pct": bacterial_grade(breadth, reads, supported, t)
           for t in bac["sensitivity_breadths"]},
    })
fields = list(records[0])
write(args.out_bacterial, fields, records)

# Viruses: high confidence when a member of the vOTU was assembled from this library.
votu_libraries = defaultdict(set)
for path in args.votus:
    for row in rows(path):
        votu_libraries[row["votu_id"]].add(row["sample_id"])
viral = []
for row in rows(args.viral):
    breadth = float(row["covered_fraction"])
    if breadth >= vir["present"]["min_breadth"]:
        grade = "high_confidence" if row["sample_id"] in votu_libraries[row["target_id"]] else "validated"
    else:
        grade = "partial" if breadth >= vir["partial"]["min_breadth"] else "none"
    viral.append({"sample_id": row["sample_id"], "target_id": row["target_id"],
                  "catalog_class": row["catalog_class"], "breadth": round(breadth, 6),
                  "reads": row["read_count"], "grade": grade})
write(args.out_viral, list(viral[0]), viral)
print(f"bacterial pairs={len(records)} viral pairs={len(viral)}")
