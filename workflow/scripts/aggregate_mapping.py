#!/usr/bin/env python3
"""Combine per-library CoverM tables into long-form coverage and library accounting tables."""
import argparse
import csv
import json
import os
import re
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--samples", required=True, nargs="+")
parser.add_argument("--results", required=True, type=Path)
parser.add_argument("--provenance", required=True, type=Path)
parser.add_argument("--manifest", required=True)
parser.add_argument("--eligibility", required=True)
parser.add_argument("--samples-table", required=True)
parser.add_argument("--bacterial", required=True, type=Path)
parser.add_argument("--viral", required=True, type=Path)
parser.add_argument("--libraries", required=True, type=Path)
args = parser.parse_args()

METRICS = {"Covered Fraction": "covered_fraction", "Read Count": "read_count", "Mean": "mean_depth",
           "Trimmed Mean": "trimmed_mean_depth", "Covered Bases": "covered_bases", "Length": "length"}


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def coverm(path):
    """CoverM columns are '<bam stem> <Metric>'; rename to plain metric names."""
    records = []
    for row in rows(path):
        target = row.get("Genome") or row.get("Contig")
        record = {"target_id": target}
        for column, value in row.items():
            for suffix, name in METRICS.items():
                if column.endswith(f" {suffix}"):
                    record[name] = value
        records.append(record)
    return records


def primary_mapped(path):
    text = Path(path).read_text()
    match = re.search(r"^(\d+) \+ \d+ mapped", text, re.M)
    return int(match.group(1)) if match else 0


manifest = {r["catalog_id"]: r for r in rows(args.manifest) if r["catalog_id"]}
eligibility = {r["sample_id"]: r for r in rows(args.eligibility)}
route = {r["sample_id"]: r["host_route"] for r in rows(args.samples_table)}


def write(path, fields, records):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    os.replace(temporary, path)


args.bacterial.parent.mkdir(parents=True, exist_ok=True)
bacterial, viral, libraries = [], [], []
for sample in args.samples:
    for record in coverm(args.results / "coverage" / "bacterial" / f"{sample}.tsv"):
        if record["target_id"] == "unmapped":
            continue
        entry = manifest.get(record["target_id"], {})
        bacterial.append({"sample_id": sample, **record, "role": entry.get("role", ""),
                          "taxonomy": entry.get("taxonomy", "")})
    for record in coverm(args.results / "coverage" / "viral" / f"{sample}.tsv"):
        if record["target_id"] != "unmapped":
            viral.append({"sample_id": sample, **record,
                          "catalog_class": "associate" if record["target_id"].startswith("VOTU")
                          else "endogenous_candidate"})
    if route[sample] == "none":
        trim = json.loads((args.provenance / "trim" / f"{sample}.json").read_text())
        source, pairs = "trimmed_regenerated", int(trim["pairs_after_fastp"])
    else:
        source, pairs = "host_depleted", int(eligibility[sample]["retained_pairs"])
    libraries.append({
        "sample_id": sample, "host_route": route[sample], "read_source": source, "input_pairs": pairs,
        "bacterial_filtered_reads": primary_mapped(args.results / "bam" / "bacterial" / f"{sample}.flagstat.txt"),
        "viral_filtered_reads": primary_mapped(args.results / "bam" / "viral" / f"{sample}.flagstat.txt"),
    })

write(args.bacterial, ["sample_id", "target_id", "role", *METRICS.values(), "taxonomy"], bacterial)
write(args.viral, ["sample_id", "target_id", "catalog_class", *METRICS.values()], viral)
write(args.libraries, ["sample_id", "host_route", "read_source", "input_pairs", "bacterial_filtered_reads",
                       "viral_filtered_reads"], libraries)
print(f"samples={len(args.samples)} bacterial_rows={len(bacterial)} viral_rows={len(viral)}")
