#!/usr/bin/env python3
"""Summarize presence-rule control behaviour from Phase-5 bacterial coverage.

Positives: every Phase-3 MAG species in the library it was assembled from (via the
catalog genome that represents it). Negatives: the 20 decoy genomes in every
library. Writes the tables used to lock the presence rule to
data/results/phase5_mapping/<version>/<scope>/presence_controls/.
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
parser.add_argument("--version", default="v2", help="catalog version")
args = parser.parse_args()
OUT = RESULTS / "phase5_mapping" / args.version / args.scope / "presence_controls"


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


manifest = rows(RESULTS / "phase4_catalog" / args.version / "bacterial_catalog.manifest.tsv")
mags = {r["species_cluster"]: r for r in rows(RESULTS / "phase3_catalog" / "mags" / "mag_catalog.tsv")
        if r["catalog_status"] == "retained" and r["mag_id"] == r["representative"]}
# Map each MAG species to the catalog genome that represents it (itself or a reference).
representative_of = {}
for entry in manifest:
    if entry["status"] == "representative":
        for member in entry["cluster_members"].split(","):
            representative_of[member] = entry["catalog_id"]
coverage = {(r["sample_id"], r["target_id"]): r for r in rows(RESULTS / "phase5_mapping" / args.version / args.scope /
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


READ_LENGTH = 150  # nominal paired-read length; trimmed reads are slightly shorter


def metrics(record):
    """Breadth, reads, read-based depth, and observed/expected breadth.

    Depth is derived from read count rather than CoverM's mean, which reports 0 for
    targets whose reads pile onto a single conserved locus.
    """
    breadth, reads = float(record["covered_fraction"]), int(record["read_count"])
    depth = reads * READ_LENGTH / int(record["length"]) if int(record["length"]) else 0.0
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

# Breadth-ratio distribution among candidate calls (>= 10% breadth), by catalog role.
calls = [(r, metrics(r)) for (s, t), r in coverage.items() if r["role"] != "decoy" and metrics(r)[0] >= 0.10]
ratio_rows = []
for low in [i / 10 for i in range(10)]:
    group = [(r, m) for r, m in calls if low <= m[3] < low + 0.1 or (low == 0.9 and m[3] >= 1)]
    ratio_rows.append([f"{low:.1f}-{low + 0.1:.1f}", len(group), sum(r["role"] == "mag" for r, _ in group)])
write("breadth_ratio_calls.tsv", ["breadth_ratio", "calls", "mag_calls"], ratio_rows)
bands = [(0, 1e-12, "0"), (1e-12, 0.01, "<1%"), (0.01, 0.05, "1-5%"), (0.05, 0.10, "5-10%"),
         (0.10, 0.20, "10-20%"), (0.20, 0.50, "20-50%"), (0.50, 1.01, ">=50%")]
real = [metrics(r)[0] for r in coverage.values() if r["role"] != "decoy"]
write("bacterial_breadth_bands.tsv", ["breadth", "genome_library_pairs"],
      [[label, sum((b == 0) if label == "0" else (lo <= b < hi) for b in real)] for lo, hi, label in bands])

# Viruses: each associate vOTU in the library its representative came from.
viral = {(r["sample_id"], r["target_id"]): r for r in rows(RESULTS / "phase5_mapping" / args.version / args.scope /
                                                           "viral_coverage.tsv")}
votus = rows(RESULTS / "phase4_catalog" / args.version / "viral_associate.manifest.tsv")
reps = {v["representative_contig"]: v["votu_id"] for v in votus}
shared = {}
for r in rows(RESULTS / "phase3_catalog" / "votus" / "associate.ani.tsv"):
    if r["qname"] != r["tname"] and r["qname"] in reps and r["tname"] in reps and float(r["pid"]) >= 95:
        shared[r["qname"]] = max(shared.get(r["qname"], 0.0), float(r["qcov"]))
source_rows = []
for v in votus:
    library = "__".join(v["representative_contig"].split("__")[:2])
    record = viral.get((library, v["votu_id"]))
    if record:
        source_rows.append([v["votu_id"], library, float(record["covered_fraction"]), record["read_count"],
                            round(shared.get(v["representative_contig"], 0.0), 1)])
write("viral_source_breadth.tsv", ["votu_id", "source_library", "breadth", "reads",
                                   "max_percent_shared_with_other_votu"],
      [r for r in sorted(source_rows, key=lambda r: r[2]) if r[2] < 0.75])
passing = sum(r[2] >= 0.75 for r in source_rows)
write("viral_source_summary.tsv", ["votus_in_scope", "source_breadth_ge_75", "below_75_sharing_ge_30pct"],
      [[len(source_rows), passing, sum(r[2] < 0.75 and r[4] >= 30 for r in source_rows)]])
