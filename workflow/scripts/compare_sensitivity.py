#!/usr/bin/env python3
"""Sensitivity analyses C1-C4 from docs/phase6_analysis_plan.md."""
import argparse
import csv
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--presence", required=True)
parser.add_argument("--subsets", required=True)
parser.add_argument("--grades", required=True)
parser.add_argument("--cohort-viral", required=True)
parser.add_argument("--nominated", required=True)
parser.add_argument("--catalog", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--sensitivity-dir", required=True, type=Path)
parser.add_argument("--out-thresholds", required=True, type=Path)
parser.add_argument("--out-host-handling", required=True, type=Path)
parser.add_argument("--out-capping", required=True, type=Path)
parser.add_argument("--out-loso", required=True, type=Path)
args = parser.parse_args()
presence = json.loads(Path(args.presence).read_text())
PRESENT = ("validated", "high_confidence")


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records):
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    os.replace(temporary, path)


def coverm(path):
    """Target -> (breadth, reads) from a CoverM table."""
    out = {}
    for r in rows(path):
        target = r.get("Genome") or r.get("Contig")
        if target == "unmapped":
            continue
        breadth = next(float(v) for k, v in r.items() if k.endswith(" Covered Fraction"))
        reads = next(int(float(v)) for k, v in r.items() if k.endswith(" Read Count"))
        out[target] = (breadth, reads)
    return out


def called(kind, breadth, reads):
    if kind == "bacterial":
        return breadth >= presence["bacteria"]["validated"]["min_breadth"] and \
            reads >= presence["bacteria"]["validated"]["min_reads"]
    return breadth >= presence["viruses"]["present"]["min_breadth"]


grades = rows(args.grades)
real = [r for r in grades if r["role"] != "decoy"]
meta = {}
with open(args.manifest, newline="") as handle:
    for r in csv.DictReader(handle):
        if r["include_primary"].lower() == "true":
            meta[f"{r['study']}__{r['library_id'].split(':')[-1]}"] = r

# C1: presence-threshold sensitivity.
thresholds = []
for column, label in (("grade_at_5pct", "5%"), ("grade", "10% (primary)"), ("grade_at_20pct", "20%")):
    calls = [r for r in real if r[column] in PRESENT]
    thresholds.append({"breadth_threshold": label, "present_pairs": len(calls),
                       "genomes": len({r["target_id"] for r in calls}),
                       "libraries": len({r["sample_id"] for r in calls}),
                       "by_study": ";".join(f"{k}:{v}" for k, v in
                                            sorted(Counter(meta[r["sample_id"]]["study"] for r in calls).items()))})
write(args.out_thresholds, thresholds)

# Reference calls from the cohort (Phase-5 read source) for comparison.
cohort = {"bacterial": defaultdict(dict), "viral": defaultdict(dict)}
for r in grades:
    cohort["bacterial"][r["sample_id"]][r["target_id"]] = (float(r["breadth"]), int(r["reads"]))
for r in rows(args.cohort_viral):
    cohort["viral"][r["sample_id"]][r["target_id"]] = (float(r["covered_fraction"]), int(r["read_count"]))
decoys = {r["target_id"] for r in grades if r["role"] == "decoy"}


def compare(sample, kind, a, b, label_a, label_b, analysis):
    targets = set(a) | set(b)
    in_a = {t for t in targets if t not in decoys and called(kind, *a.get(t, (0, 0)))}
    in_b = {t for t in targets if t not in decoys and called(kind, *b.get(t, (0, 0)))}
    union = in_a | in_b
    return {"analysis": analysis, "sample_id": sample, "domain": kind, "comparison": f"{label_a} vs {label_b}",
            f"present_{label_a}": len(in_a), f"present_{label_b}": len(in_b), "both": len(in_a & in_b),
            f"only_{label_a}": len(in_a - in_b), f"only_{label_b}": len(in_b - in_a),
            "jaccard": round(len(in_a & in_b) / len(union), 4) if union else 1.0}


subsets = rows(args.subsets)
host_rows, cap_rows = [], []
for s in subsets:
    sample = s["sample_id"]
    for kind in ("bacterial", "viral"):
        path = lambda variant: args.sensitivity_dir / variant / "coverage" / kind / f"{sample}.tsv"
        if s["analysis"] == "host_handling":
            host_rows.append(compare(sample, kind, cohort[kind][sample], coverm(path("full_reads")),
                                     "host_depleted", "full_reads", "host_handling"))
        else:
            capped = cohort[kind][sample] if s["stratum"] == "reference_free" else coverm(path("full_reads"))
            cap_rows.append(compare(sample, kind, capped, coverm(path("uncapped")), "capped", "uncapped", "capping"))
write(args.out_host_handling, host_rows)
write(args.out_capping, cap_rows)

# C4: presences whose catalog genome was nominated only by one study's libraries.
catalog_of = {}
with open(args.catalog, newline="") as handle:
    for e in csv.DictReader(handle, delimiter="\t"):
        if e["status"] == "representative":
            for member in e["cluster_members"].split(","):
                catalog_of[member] = e["catalog_id"]
nominating = defaultdict(set)
for r in rows(args.nominated):
    target = catalog_of.get(r["genome_id"])
    if target:
        nominating[target] |= {meta[s]["study"] for s in r["library_list"].split(",") if s in meta}
present = [r for r in real if r["grade"] in PRESENT]
mag_only = {r["target_id"] for r in real} - set(nominating)
loso = []
for study in sorted({m["study"] for m in meta.values()}):
    single = {t for t, studies in nominating.items() if studies == {study}}
    dependent = [r for r in present if r["target_id"] in single]
    loso.append({"study": study, "genomes_nominated_only_by_study": len(single),
                 "dependent_presences": len(dependent),
                 "dependent_presences_in_other_studies": sum(meta[r["sample_id"]]["study"] != study for r in dependent),
                 "share_of_all_presences": round(len(dependent) / len(present), 4) if present else 0})
loso.append({"study": "MAG_only_genomes (no sylph nomination)", "genomes_nominated_only_by_study": len(mag_only),
             "dependent_presences": sum(r["target_id"] in mag_only for r in present),
             "dependent_presences_in_other_studies": "",
             "share_of_all_presences": round(sum(r["target_id"] in mag_only for r in present) / len(present), 4)})
write(args.out_loso, loso)
print(f"thresholds={len(thresholds)} host_handling={len(host_rows)} capping={len(cap_rows)} loso={len(loso)}")
