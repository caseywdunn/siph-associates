#!/usr/bin/env python3
"""Calibrate the eukaryotic read-level presence rule against negative-control lineages.

Reads data/results/eukaryote_gate/read_lineages.tsv (competitive SSU remap, LCA
assignment) and the assembled-SSU classifications, and writes calibration tables
to data/results/eukaryote_gate/calibration/.
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from eukaryote_units import unit  # noqa: E402

GATE = ROOT / "data" / "results" / "eukaryote_gate"
OUT = GATE / "calibration"
BANDS = (97, 99)
THRESHOLDS = (2, 5, 10, 25, 50)


def write(name, header, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n# {name}")
    for row in [header, *rows[:40]]:
        print("\t".join(str(v) for v in row))


pairs = defaultdict(Counter)  # (sample, unit) -> band -> cumulative read pairs
kinds = {}
with open(GATE / "read_lineages.tsv", newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        u = unit(row["lineage"])
        if u is None:
            continue
        kinds[u[1]] = u[0]
        for band in BANDS:
            if int(row["min_identity_band"]) >= band:
                pairs[(row["sample_id"], u[1])][band] += int(row["read_pairs"])

# Negative controls: background under identical filters.
rows = []
for band in BANDS:
    for control in ("Embryophyta", "Insecta", "Mammalia"):
        values = [c[band] for (s, u), c in pairs.items() if u == control]
        rows.append([f">={band}%", control, sum(v > 0 for v in values), max(values, default=0)])
write("negative_controls.tsv", ["identity", "control", "libraries_with_pairs", "max_pairs_in_a_library"], rows)
background = {band: max([c[band] for (s, u), c in pairs.items() if kinds.get(u) == "negative_control"], default=0)
              for band in BANDS}

# Positives: non-host clades assembled as full-length SSUs in the same library.
positives = []
for classification in sorted((GATE / "ssu").glob("*.ssu_classification.csv")):
    sample = classification.name.removesuffix(".ssu_classification.csv")
    with classification.open(newline="") as handle:
        for row in csv.DictReader(handle):
            u = unit(row["taxonomy"])
            if u and u[0] in ("metazoan", "non_metazoan"):
                positives.append((sample, u[1], float(row["%id"] or 0)))
positive_rows = []
for sample, u, identity in sorted(set(positives)):
    c = pairs.get((sample, u), Counter())
    positive_rows.append([sample, u, identity, c[97], c[99]])
write("assembled_positives.tsv", ["sample_id", "unit", "assembled_identity", "pairs_ge97", "pairs_ge99"],
      positive_rows)

# Candidate detections by threshold, excluding host, cnidarians, human, and negative controls.
grid = []
for band in BANDS:
    for n in THRESHOLDS:
        detected = [(s, u) for (s, u), c in pairs.items() if kinds.get(u) in ("metazoan", "non_metazoan")
                    and c[band] >= n]
        recovered = sum(pairs.get((s, u), Counter())[band] >= n for s, u, _ in set(positives))
        negatives = sum(c[band] >= n for (s, u), c in pairs.items() if kinds.get(u) == "negative_control")
        grid.append([f">={band}%", n, len(detected), len({s for s, _ in detected}), negatives,
                     f"{recovered}/{len(set(positives))}"])
write("threshold_grid.tsv", ["identity", "min_pairs", "detections", "libraries", "negative_control_detections",
                             "assembled_positives_recovered"], grid)

for band in BANDS:
    n = background[band] + 1
    detected = [(s, u) for (s, u), c in pairs.items() if kinds.get(u) in ("metazoan", "non_metazoan")
                and c[band] >= n]
    write(f"detections_above_background_ge{band}.tsv", ["unit", "kind", "libraries"],
          [[u, kinds[u], k] for u, k in Counter(u for _, u in detected).most_common()])
cnidarian = [c[99] for (s, u), c in pairs.items() if u == "Cnidaria"]
print(f"\nbackground (max negative-control pairs): {background}; non-host cnidarian pairs >=99%: "
      f"libraries {sum(v > 0 for v in cnidarian)}, max {max(cnidarian, default=0)}")
