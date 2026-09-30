#!/usr/bin/env python3
"""Freeze the Phase-6 sensitivity subsets from library metadata only (never results).

Host handling (C2): 5 P. physalis and 5 N. septata assembly-eligible libraries,
one seeded draw per input-depth quintile within each species, remapped with full
trimmed reads. Capping (C3): 8 libraries above the 200 M-pair cap, 4 reference-free
and 4 reference-bearing, seeded; each is mapped capped and uncapped with full
trimmed reads so the comparison isolates capping.
"""
import csv
import json
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
settings = json.loads((ROOT / "config" / "phase6_analysis.json").read_text())["sensitivity"]
cap = 200_000_000
with open(ROOT / "config" / "phase3_assembly.tsv", newline="") as handle:
    members = {r["sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
with open(ROOT / "config" / "samples.tsv", newline="") as handle:
    samples = {r["sample_id"]: r for r in csv.DictReader(handle, delimiter="\t")
               if r["include_primary"].lower() == "true"}

rows = []
host = settings["host_handling_subset"]
generator = random.Random(host["seed"])
for route in ("P_physalis", "N_septata"):
    pool = sorted((s for s, m in members.items() if m["host_route"] == route and m["assembly_eligible"] == "true"),
                  key=lambda s: (int(members[s]["retained_pairs"]), s))
    n = host[route]
    for q in range(n):
        stratum = pool[q * len(pool) // n:(q + 1) * len(pool) // n]
        rows.append({"sample_id": generator.choice(stratum), "analysis": "host_handling",
                     "variants": "full_reads", "stratum": f"{route}:depth_quintile_{q + 1}"})

capping = settings["capping_subset"]
generator = random.Random(capping["seed"])
for route_group, n in (("reference_free", capping["libraries"] // 2), ("reference_bearing", capping["libraries"] // 2)):
    pool = sorted(s for s, r in samples.items() if int(r["read_pairs"]) > cap
                  and (r["host_route"] == "none") == (route_group == "reference_free"))
    for sample in sorted(generator.sample(pool, n)):
        variants = "uncapped" if route_group == "reference_free" else "full_reads;uncapped"
        rows.append({"sample_id": sample, "analysis": "capping", "variants": variants, "stratum": route_group})

output = ROOT / "config" / "phase6_subsets.tsv"
with output.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["sample_id", "analysis", "variants", "stratum"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(output.read_text())
