#!/usr/bin/env python3
"""Profile eukaryotic SSU evidence for the eukaryote gate (evidence review only).

Reads the per-library extractions from extract_eukaryote_evidence.py and writes
summary tables to data/results/eukaryote_gate/summary/.
"""
from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
EVIDENCE = ROOT / "data" / "results" / "eukaryote_gate" / "evidence"
OUT = ROOT / "data" / "results" / "eukaryote_gate" / "summary"


def category(lineage):
    ranks = lineage.split(";")
    if ranks[0] != "Eukaryota":
        return "prokaryote_or_organelle", ranks[1] if len(ranks) > 1 else ""
    if "Homo sapiens" in lineage or ";Homo;" in lineage + ";":
        return "human", "Homo"
    if "Metazoa" in ranks:
        after = ranks[ranks.index("Metazoa") + 1:]
        phylum = next((r for r in after if r not in ("Animalia", "BCP clade", "Bilateria", "Protostomia",
                                                    "Deuterostomia", "Spiralia", "Ecdysozoa")), "Metazoa")
        if "Siphonophorae" in ranks:
            return "host_siphonophore", "Siphonophorae"
        if phylum == "Cnidaria":
            return "other_cnidarian", ranks[ranks.index("Cnidaria") + 1] if "Cnidaria" in ranks[:-1] else "Cnidaria"
        return "non_cnidarian_metazoan", phylum
    major = ranks[1] if len(ranks) > 1 else "Eukaryota"
    return "non_metazoan_eukaryote", major


def write(name, header, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / name, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n# {name}")
    for row in [header, *rows[:40]]:
        print("\t".join(str(v) for v in row))


ssus, missing = [], []
for directory in sorted(EVIDENCE.iterdir()):
    table = directory / f"{directory.name}.phyloFlash.extractedSSUclassifications.csv"
    if not table.is_file():
        missing.append(directory.name)
        continue
    with table.open(newline="") as handle:
        for row in csv.DictReader(handle):
            group, clade = category(row["taxonomy"])
            ssus.append({"sample": directory.name, "otu": row["OTU"], "group": group, "clade": clade,
                         "lineage": row["taxonomy"], "identity": float(row["%id"] or 0),
                         "length": int(row["alnlen"] or 0), "read_cov": int(float(row["read_cov"] or 0))})
libraries = sorted(d.name for d in EVIDENCE.iterdir())
print(f"libraries with evidence: {len(libraries)}; missing classification tables: {len(missing)}")

groups = defaultdict(list)
for s in ssus:
    groups[s["group"]].append(s)
write("assembled_ssu_groups.tsv", ["group", "assembled_ssus", "libraries", "median_identity", "median_read_cov"],
      [[g, len(v), len({s["sample"] for s in v}), sorted(x["identity"] for x in v)[len(v) // 2],
        sorted(x["read_cov"] for x in v)[len(v) // 2]] for g, v in sorted(groups.items())])

clades = defaultdict(list)
for s in ssus:
    if s["group"] in ("non_cnidarian_metazoan", "other_cnidarian", "non_metazoan_eukaryote", "human"):
        clades[(s["group"], s["clade"])].append(s)
write("eukaryote_clades.tsv", ["group", "clade", "assembled_ssus", "libraries", "identity_min", "identity_median",
                               "example_lineage"],
      [[g, c, len(v), len({s["sample"] for s in v}), min(x["identity"] for x in v),
        sorted(x["identity"] for x in v)[len(v) // 2], v[0]["lineage"].split(";")[-3:]]
       for (g, c), v in sorted(clades.items(), key=lambda kv: -len({s["sample"] for s in kv[1]}))])

# Identity bands: how many non-host SSUs could be named at genus/species level against SILVA.
bands = Counter()
for s in ssus:
    if s["group"] in ("non_cnidarian_metazoan", "non_metazoan_eukaryote", "other_cnidarian"):
        bands["<90%" if s["identity"] < 90 else "90-97%" if s["identity"] < 97 else ">=97%"] += 1
write("nonhost_identity_bands.tsv", ["identity_to_silva", "assembled_ssus"],
      [[b, bands[b]] for b in ("<90%", "90-97%", ">=97%")])

# Read-level signal outside Metazoa (NTU table is truncated at seven ranks).
reads = defaultdict(Counter)
for directory in sorted(EVIDENCE.iterdir()):
    table = directory / f"{directory.name}.phyloFlash.NTUfull_abundance.csv"
    if not table.is_file():
        continue
    for line in table.read_text().splitlines():
        lineage, _, count = line.rpartition(",")
        if lineage.startswith("Eukaryota"):
            ranks = lineage.split(";")
            key = "Metazoa (unresolved)" if "Metazoa" in ranks else ("unresolved" if ranks[1].startswith("(")
                                                                     else ranks[1])
            reads[key][directory.name] += int(count)
write("read_level_eukaryotes.tsv", ["group", "total_reads", "libraries_with_reads", "libraries_with_ge_100_reads"],
      [[k, sum(v.values()), len(v), sum(n >= 100 for n in v.values())]
       for k, v in sorted(reads.items(), key=lambda kv: -sum(kv[1].values()))])
