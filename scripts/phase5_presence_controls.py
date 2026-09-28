#!/usr/bin/env python3
"""Summarize presence-rule control behaviour from Phase-5 bacterial coverage.

Positives: every Phase-3 MAG species in the library it was assembled from (via the
catalog genome that represents it). Negatives: the 20 decoy genomes in every
library. Writes the tables used to lock the presence rule to
data/results/phase5_mapping/<scope>/presence_controls/.
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "data" / "results"
parser = argparse.ArgumentParser()
parser.add_argument("--scope", choices=("pilot", "cohort"), default="cohort")
args = parser.parse_args()
OUT = RESULTS / "phase5_mapping" / args.scope / "presence_controls"


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name, header, records):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(records)
    print(f"\n# {name}")
    for record in [header, *records]:
        print("\t".join(str(v) for v in record))


manifest = rows(RESULTS / "phase4_catalog" / "v1" / "bacterial_catalog.manifest.tsv")
mags = {r["species_cluster"]: r for r in rows(RESULTS / "phase3_catalog" / "mags" / "mag_catalog.tsv")
        if r["catalog_status"] == "retained" and r["mag_id"] == r["representative"]}
# Map each MAG species to the catalog genome that represents it (itself or a reference).
representative_of = {}
for entry in manifest:
    if entry["status"] == "representative":
        for member in entry["cluster_members"].split(","):
            representative_of[member] = entry["catalog_id"]
coverage = {(r["sample_id"], r["target_id"]): r for r in rows(RESULTS / "phase5_mapping" / args.scope /
                                                             "bacterial_coverage.tsv")}
samples = sorted({s for s, _ in coverage})
decoys = sorted(e["catalog_id"] for e in manifest if e["role"] == "decoy" and e["catalog_id"])

positives = []
for species, mag in sorted(mags.items()):
    target = representative_of.get(species)
    key = (mag["sample_id"], target)
    if key in coverage:
        positives.append((species, *key, coverage[key]))
negatives = [(s, d, coverage[(s, d)]) for s in samples for d in decoys if (s, d) in coverage]


def expected_breadth(depth):
    """Breadth expected for uniform coverage at this mean depth (Poisson, 1 - e^-depth)."""
    return 1 - math.exp(-depth)


def metrics(record):
    breadth, reads, depth = float(record["covered_fraction"]), int(record["read_count"]), float(record["mean_depth"])
    expected = expected_breadth(depth)
    return breadth, reads, depth, (breadth / expected if expected > 0 else 0.0)


write("positives.tsv", ["mag_species", "source_library", "catalog_id", "breadth", "reads", "mean_depth",
                        "breadth_ratio"],
      [[sp, s, t, *[round(v, 4) for v in metrics(r)]] for sp, s, t, r in positives])
decoy_hits = [(s, d, *metrics(r)) for s, d, r in negatives if int(r["read_count"]) > 0]
write("decoy_recruitment.tsv", ["sample_id", "decoy_id", "breadth", "reads", "mean_depth", "breadth_ratio"],
      [[s, d, *[round(v, 4) for v in m]] for s, d, *m in sorted(decoy_hits, key=lambda x: -x[2])[:50]])

grid = []
for min_breadth in (0.01, 0.05, 0.10, 0.20, 0.50):
    for min_reads in (10, 50, 100):
        for min_ratio in (0.0, 0.5):
            passes = lambda r: (metrics(r)[0] >= min_breadth and metrics(r)[1] >= min_reads
                                and metrics(r)[3] >= min_ratio)
            tp = sum(passes(r) for *_, r in positives)
            fp = sum(passes(r) for *_, r in negatives)
            grid.append([min_breadth, min_reads, min_ratio, f"{tp}/{len(positives)}",
                         f"{fp}/{len(negatives)}", len({s for s, _, r in negatives if passes(r)})])
write("rule_grid.tsv", ["min_breadth", "min_reads", "min_breadth_ratio", "positives_passing",
                        "decoy_library_pairs_passing", "libraries_with_decoy_pass"], grid)
