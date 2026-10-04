#!/usr/bin/env python3
"""Primary analyses B1, B4, B5: incidence and grade summaries by host and study; CRISPR-supported links."""
import argparse
import csv
import os
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--bacterial", required=True)
parser.add_argument("--grade-column", default="grade", choices=("grade", "grade_at_5pct", "grade_at_20pct"))
parser.add_argument("--viral", required=True)
parser.add_argument("--contamination", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--samples", required=True)
parser.add_argument("--links", required=True)
parser.add_argument("--out-bacterial", required=True, type=Path)
parser.add_argument("--out-viral", required=True, type=Path)
parser.add_argument("--out-grades", required=True, type=Path)
parser.add_argument("--out-links", required=True, type=Path)
args = parser.parse_args()
PRESENT = ("validated", "high_confidence")


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


meta = {}
with open(args.manifest, newline="") as handle:
    for r in csv.DictReader(handle):
        if r["include_primary"].lower() == "true":
            meta[f"{r['study']}__{r['library_id'].split(':')[-1]}"] = r
route = {r["sample_id"]: ("reference_free" if r["host_route"] == "none" else "reference_bearing")
         for r in rows(args.samples)}
contamination = {r["target_id"]: r for r in rows(args.contamination)}


def incidence(table, kind):
    by_target = defaultdict(list)
    for r in table:
        by_target[r["target_id"]].append(r)
    out = []
    for target, group in sorted(by_target.items()):
        present = [r for r in group if r["grade"] in PRESENT]
        species = Counter(meta[r["sample_id"]]["species_current"] for r in present)
        studies = Counter(meta[r["sample_id"]]["study"] for r in present)
        record = {"target_id": target,
                  "present_libraries": len(present),
                  "high_confidence_libraries": sum(r["grade"] == "high_confidence" for r in group),
                  "host_species": len(species),
                  "by_host_species": ";".join(f"{k}:{v}" for k, v in species.most_common()),
                  "by_study": ";".join(f"{k}:{v}" for k, v in sorted(studies.items()))}
        if kind == "bacterial":
            c = contamination.get(target, {})
            record.update({"role": group[0]["role"], "taxonomy": group[0]["taxonomy"],
                           "trace_libraries": sum(r["grade"] == "trace" for r in group),
                           "divergent_strain_calls": sum(r["divergent_strain"] == "true" for r in present),
                           "contamination_status": c.get("status", "below_minimum_presences"),
                           "probable_contaminant": c.get("probable_contaminant", "false"),
                           "contamination_q": c.get("q_value", ""),
                           "identity_annotation": c.get("identity_annotation", "")})
        else:
            record.update({"catalog_class": group[0]["catalog_class"],
                           "partial_libraries": sum(r["grade"] == "partial" for r in group)})
        out.append(record)
    return out


bacterial = rows(args.bacterial)
# Use the selected threshold consistently for counts, confidence, and linkage.
for row in bacterial:
    row["grade"] = row[args.grade_column]
    if args.grade_column != "grade":
        row["divergent_strain"] = str(
            row["grade"] in PRESENT and float(row["breadth_ratio"]) < 0.5
        ).lower()
viral = rows(args.viral)
b = incidence(bacterial, "bacterial")
v = incidence(viral, "viral")
write(args.out_bacterial, list(b[0]), b)
write(args.out_viral, list(v[0]), v)

# Grade counts by host route and domain (reference-free assembly used a subsample).
grades = []
for kind, table in (("bacterial", [r for r in bacterial if r["role"] != "decoy"]), ("viral", viral)):
    counts = Counter((route[r["sample_id"]], r["grade"]) for r in table)
    for (r_, g), n in sorted(counts.items()):
        grades.append({"domain": kind, "host_route": r_, "grade": g, "pairs": n})
write(args.out_grades, ["domain", "host_route", "grade", "pairs"], grades)

# B5: phage-host links only where a CRISPR spacer supports them; co-occurrence is descriptive.
present_b = defaultdict(set)
for r in bacterial:
    if r["grade"] in PRESENT:
        present_b[r["target_id"]].add(r["sample_id"])
present_v = defaultdict(set)
for r in viral:
    if r["grade"] in PRESENT:
        present_v[r["target_id"]].add(r["sample_id"])
links = []
for r in rows(args.links):
    links.append({**r, "votu_present_libraries": len(present_v[r["votu_id"]]),
                  "host_present_libraries": len(present_b[r["host_catalog_id"]]),
                  "co_present_libraries": len(present_v[r["votu_id"]] & present_b[r["host_catalog_id"]])})
write(args.out_links, list(links[0]) if links else ["votu_id"], links)
print(f"bacterial targets={len(b)} viral targets={len(v)} links={len(links)}")
