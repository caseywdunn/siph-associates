#!/usr/bin/env python3
"""Tabulate per-library assembly-branch products, including documented exclusions."""
import argparse
import csv
import os
import sys
import tarfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import fasta_lengths, fasta_summary, sha256

FIELDS = [
    "sample_id", "host_route", "strategy", "retained_pairs", "assembly_status",
    "contigs", "assembly_bases", "n50", "max_contig", "assembly_sha256",
    "rrna_markers", "bacterial_rrna", "archaeal_rrna", "eukaryotic_rrna",
    "viral_contigs", "viral_bases", "checkv_complete", "checkv_high_quality", "checkv_medium_quality",
    "binnable_contigs", "raw_bins", "binned_bases", "checkm2_status", "assessed_bins", "unassessed_bins",
    "medium_quality_or_better_bins", "high_quality_bins", "assembly_seconds", "assembly_max_rss_mb",
]

parser = argparse.ArgumentParser()
parser.add_argument("--membership", required=True)
parser.add_argument("--results", required=True, type=Path)
parser.add_argument("--benchmarks", required=True, type=Path)
parser.add_argument("--metabat-min-contig", required=True, type=int)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize(sample):
    path = lambda kind, suffix: args.results / kind / f"{sample}{suffix}"
    contigs = fasta_summary(path("assemblies", ".fasta"))
    lengths = fasta_lengths(path("assemblies", ".fasta"))
    kingdoms = Counter()
    for line in path("markers", ".gff").read_text().splitlines():
        kingdoms[line.rsplit("kingdom=", 1)[-1]] += 1
    viruses = fasta_lengths(path("viruses", ".fna"))
    tiers = Counter(row["checkv_quality"] for row in rows(path("checkv", ".tsv")))
    bin_bases, bins = 0, 0
    with tarfile.open(path("bins", ".tar.gz"), "r:gz") as archive:
        for member in archive.getmembers():
            if member.name.endswith(".fa"):
                bins += 1
                bin_bases += sum(len(line.strip()) for line in archive.extractfile(member)
                                 if not line.startswith(b">"))
    quality = rows(path("checkm2", ".tsv"))
    timing = rows(args.benchmarks / f"{sample}.tsv")[-1]
    return {
        "assembly_status": "assembled" if contigs["contigs"] else "no_contigs",
        "contigs": contigs["contigs"], "assembly_bases": contigs["total_bases"],
        "n50": contigs["n50"], "max_contig": contigs["max_contig"],
        "assembly_sha256": sha256(path("assemblies", ".fasta")),
        "rrna_markers": sum(kingdoms.values()), "bacterial_rrna": kingdoms["bac"],
        "archaeal_rrna": kingdoms["arc"], "eukaryotic_rrna": kingdoms["euk"],
        "viral_contigs": len(viruses), "viral_bases": sum(viruses),
        "checkv_complete": tiers["Complete"], "checkv_high_quality": tiers["High-quality"],
        "checkv_medium_quality": tiers["Medium-quality"],
        "binnable_contigs": sum(length >= args.metabat_min_contig for length in lengths),
        "raw_bins": bins, "binned_bases": bin_bases,
        "checkm2_status": path("checkm2", ".status").read_text().strip(),
        "assessed_bins": len(quality), "unassessed_bins": bins - len(quality),
        "medium_quality_or_better_bins": sum(
            float(row["Completeness"]) >= 50 and float(row["Contamination"]) < 10 for row in quality),
        "high_quality_bins": sum(
            float(row["Completeness"]) >= 90 and float(row["Contamination"]) < 5 for row in quality),
        "assembly_seconds": timing["s"], "assembly_max_rss_mb": timing["max_rss"],
    }


table = []
for member in rows(args.membership):
    row = dict.fromkeys(FIELDS, "")
    row.update({key: member[key] for key in ("sample_id", "host_route", "strategy", "retained_pairs")})
    if member["assembly_eligible"] == "true":
        row.update(summarize(member["sample_id"]))
    else:
        row["assembly_status"] = "excluded_below_floor"
    table.append(row)

temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(table)
os.replace(temporary, args.output)
print("".join(f"{status}={count}\n" for status, count in Counter(r["assembly_status"] for r in table).items()), end="")
