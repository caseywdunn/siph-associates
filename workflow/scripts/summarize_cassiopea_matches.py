#!/usr/bin/env python3
"""Summarize Cassiopea 16S matches to the siphonophore Mycoplasmatales MAG 16S genes.

Each denoised variant (ZOTU) or PacBio read is assigned to its best-matching MAG species and an
identity band: >= species threshold, >= genus-level threshold, or below (down to the reporting minimum).
"""
import argparse
import csv
import json
import os
import statistics
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--runs", required=True)
parser.add_argument("--table", required=True)
parser.add_argument("--zotus", required=True)
parser.add_argument("--pacbio", required=True)
parser.add_argument("--pacbio-reads", required=True)
parser.add_argument("--samples", required=True, type=Path)
parser.add_argument("--types", required=True, type=Path)
parser.add_argument("--pacbio-summary", required=True, type=Path)
args = parser.parse_args()
settings = json.loads(Path(args.config).read_text())["match"]
species_min, genus_min = settings["bands"]["same_species_full_length"], settings["bands"]["genus_level"]


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    os.replace(temporary, path)


def band(identity):
    if identity >= species_min:
        return f">={species_min:g}"
    if identity >= genus_min:
        return f"{genus_min:g}-{species_min:g}"
    return f"{settings['min_identity']:g}-{genus_min:g}"


runs = {r["run"]: r for r in rows(args.runs)}
meta_fields = ["site", "sample_type", "preservation", "animal", "sample_alias"]

# Illumina: reads per sample in ZOTUs matching each MAG species, by identity band.
match = {r["zotu"].split(";")[0]: (r["mag_16s"].split("|")[0], float(r["identity"])) for r in rows(args.zotus)}
with open(args.table) as handle:
    header = handle.readline().rstrip("\n").split("\t")
    table = [line.rstrip("\n").split("\t") for line in handle]
samples = header[1:]
totals = {s: sum(int(r[i + 1]) for r in table) for i, s in enumerate(samples)}
reads = defaultdict(int)
best = {}
for r in table:
    zotu = r[0]
    if zotu not in match:
        continue
    species, identity = match[zotu]
    for i, s in enumerate(samples):
        if int(r[i + 1]):
            reads[(s, species, band(identity))] += int(r[i + 1])
            best[(s, species)] = max(best.get((s, species), 0), identity)
records = []
for (s, species, b), n in sorted(reads.items()):
    records.append({"run": s, **{k: runs[s][k] for k in meta_fields}, "species_cluster": species, "band": b,
                    "reads": n, "sample_reads": totals[s], "relative_abundance": round(n / totals[s], 6),
                    "best_identity": best[(s, species)]})
write(args.samples, records, ["run"] + meta_fields + ["species_cluster", "band", "reads", "sample_reads",
                                                      "relative_abundance", "best_identity"])

# By sample type: how many samples carry each species and band, and at what abundance.
by_type = defaultdict(list)
for r in records:
    by_type[(runs[r["run"]]["sample_type"], r["species_cluster"], r["band"])].append(r["relative_abundance"])
type_counts = defaultdict(int)
for s in samples:
    type_counts[runs[s]["sample_type"]] += 1
summary = [{"sample_type": t, "species_cluster": sp, "band": b, "samples_with_match": len(v),
            "samples_of_type": type_counts[t], "median_relative_abundance": round(statistics.median(v), 6),
            "max_relative_abundance": round(max(v), 6)}
           for (t, sp, b), v in sorted(by_type.items())]
write(args.types, summary, ["sample_type", "species_cluster", "band", "samples_with_match", "samples_of_type",
                            "median_relative_abundance", "max_relative_abundance"])

# PacBio full-length reads: reads per run matching each species, by band.
in_range = {r["run"]: int(r["reads_in_length_range"]) for r in rows(args.pacbio_reads)}
pacbio = defaultdict(list)
for r in rows(args.pacbio):
    pacbio[(r["run"], r["mag_16s"].split("|")[0], band(float(r["identity"])))].append(float(r["identity"]))
pacbio_records = [{"run": run, **{k: runs[run][k] for k in meta_fields}, "species_cluster": sp, "band": b,
                   "reads": len(v), "reads_in_length_range": in_range[run],
                   "fraction": round(len(v) / in_range[run], 6), "best_identity": max(v)}
                  for (run, sp, b), v in sorted(pacbio.items())]
for run in in_range:
    if not any(r["run"] == run for r in pacbio_records):
        pacbio_records.append({"run": run, **{k: runs[run][k] for k in meta_fields}, "species_cluster": "none",
                               "band": "", "reads": 0, "reads_in_length_range": in_range[run], "fraction": 0,
                               "best_identity": ""})
write(args.pacbio_summary, pacbio_records, ["run"] + meta_fields + ["species_cluster", "band", "reads",
                                                                   "reads_in_length_range", "fraction",
                                                                   "best_identity"])
print(f"illumina samples={len(samples)} zotus_matched={len(match)} sample_matches={len(records)} "
      f"pacbio_runs={len(in_range)} pacbio_matches={sum(len(v) for v in pacbio.values())}")
