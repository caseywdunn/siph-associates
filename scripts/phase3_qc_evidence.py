#!/usr/bin/env python3
"""Tabulate the evidence behind the Phase-3 MAG and viral QC decisions.

prepare: write viral query sets for host-genome alignment and extract
         medium-quality-or-better bins for ANI clustering.
summarize: read accepted per-library outputs plus the minimap2 and skani
           results and write the tables cited in docs/phase3_qc_decisions.md.
Run both through scripts/phase3_qc_evidence.sbatch.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
COHORT = ROOT / "data" / "results" / "phase3_cohort"
OUT = COHORT / "qc_evidence"
QUERY_SETS = ("nanomia_free", "other_free", "septata_depleted")


def table(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def library(name):
    return "__".join(name.split("__")[:2])


def species():
    with open(ROOT / "manifest.csv", newline="") as handle:
        return {f"{row['study']}__{row['library_id'].split(':')[-1]}": row["species_current"]
                for row in csv.DictReader(handle)}


def routes():
    return {row["sample_id"]: row["host_route"] for row in table(ROOT / "config" / "samples.tsv")}


def medium_or_better(row):
    return float(row["Completeness"]) >= 50 and float(row["Contamination"]) < 10


def prepare():
    OUT.mkdir(exist_ok=True)
    names, host = species(), routes()
    handles = {key: open(OUT / f"{key}.fna", "w") for key in QUERY_SETS}
    for path in sorted(glob.glob(str(COHORT / "viruses" / "*.fna"))):
        sample = Path(path).stem
        if host[sample] == "N_septata":
            key = "septata_depleted"
        elif host[sample] == "none":
            key = "nanomia_free" if names[sample].startswith("Nanomia") else "other_free"
        else:
            continue
        handles[key].write(Path(path).read_text())
    for handle in handles.values():
        handle.close()
    keep = {row["Name"] for path in glob.glob(str(COHORT / "checkm2" / "*.tsv"))
            for row in table(path) if medium_or_better(row)}
    bins = OUT / "medium_bins"
    bins.mkdir(exist_ok=True)
    for path in sorted(glob.glob(str(COHORT / "bins" / "*.tar.gz"))):
        with tarfile.open(path) as archive:
            for member in archive.getmembers():
                if member.name[:-3] in keep:
                    (bins / member.name).write_bytes(archive.extractfile(member).read())
    print(f"medium-or-better bins: {len(keep)}")
    # One representative of the recurrent 5,000-bp DTR element, for the P. physalis check.
    element = sorted(row["contig_id"] for path in glob.glob(str(COHORT / "checkv" / "*.tsv"))
                     for row in table(path) if row["contig_length"] == "5000" and row["checkv_quality"] == "Complete")
    with open(OUT / "physalia_5kb.fa", "w") as handle:
        keep_record = False
        for line in open(COHORT / "viruses" / f"{library(element[0])}.fna"):
            if line.startswith(">"):
                keep_record = line[1:].split()[0] == element[0]
            if keep_record:
                handle.write(line)


def write(name, header, rows):
    with open(OUT / name, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"\n# {name}")
    for row in [header, *rows]:
        print("\t".join(str(value) for value in row))


def union_find(names, edges):
    parent = {name: name for name in names}

    def find(name):
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    for a, b in edges:
        parent[find(a)] = find(b)
    groups = defaultdict(list)
    for name in names:
        groups[find(name)].append(name)
    return sorted(groups.values(), key=len, reverse=True)


def aligned_fraction(paf, lengths):
    intervals = defaultdict(list)
    for line in open(paf):
        fields = line.split("\t")
        intervals[fields[0]].append((int(fields[2]), int(fields[3])))
    fractions = {}
    for query, spans in intervals.items():
        spans.sort()
        total, (start, end) = 0, spans[0]
        for s, e in spans[1:]:
            if s > end:
                total, start, end = total + end - start, s, e
            else:
                end = max(end, e)
        fractions[query] = (total + end - start) / lengths[query]
    return fractions


def summarize():
    strategy = {row["sample_id"]: row["strategy"] for row in table(COHORT / "assembly_summary.tsv")}

    # MAG completeness x contamination grid.
    bins = [row for path in glob.glob(str(COHORT / "checkm2" / "*.tsv")) for row in table(path)]
    contamination = (5, 10, 20, None)
    write("mag_quality_grid.tsv", ["completeness_min", *[f"contamination_lt_{c or 'any'}" for c in contamination]],
          [[c, *[sum(float(b["Completeness"]) >= c and (x is None or float(b["Contamination"]) < x)
                     for b in bins) for x in contamination]] for c in (90, 70, 50, 30, 0)])
    medium = [b for b in bins if medium_or_better(b)]
    write("mag_medium_by_strategy.tsv", ["strategy", "bins", "medium_or_better", "libraries_with_medium"],
          [[s, sum(strategy[library(b["Name"])] == s for b in bins),
            sum(strategy[library(b["Name"])] == s for b in medium),
            len({library(b["Name"]) for b in medium if strategy[library(b["Name"])] == s})]
           for s in ("reference_depleted", "fixed_effort_trimmed")])

    # MAG species clusters under alternative ANI / alignment-fraction thresholds.
    ani = table(OUT / "medium_bins_skani.tsv")
    names = sorted(path.name[:-3] for path in (OUT / "medium_bins").glob("*.fa"))
    stem = lambda path: Path(path).name[:-3]
    rows = []
    for threshold, fraction in ((95, 10), (95, 30), (95, 50), (97, 50), (99, 50)):
        edges = [(stem(r["Ref_file"]), stem(r["Query_file"])) for r in ani
                 if float(r["ANI"]) >= threshold
                 and min(float(r["Align_fraction_ref"]), float(r["Align_fraction_query"])) >= fraction]
        groups = union_find(names, edges)
        rows.append([threshold, fraction, len(groups), ",".join(str(len(g)) for g in groups[:5])])
    write("mag_species_clusters.tsv", ["ani_min", "align_fraction_min", "clusters", "largest_clusters"], rows)
    bands = Counter("<90" if a < 90 else "90-95" if a < 95 else "95-97" if a < 97 else "97-99" if a < 99 else ">=99"
                    for a in (float(r["ANI"]) for r in ani))
    write("mag_ani_bands.tsv", ["ani_band", "pairs"], [[b, bands[b]] for b in ("<90", "90-95", "95-97", "97-99", ">=99")])

    # Viral contigs: geNomad evidence joined to CheckV quality.
    checkv = {row["contig_id"]: row for path in glob.glob(str(COHORT / "checkv" / "*.tsv")) for row in table(path)}
    viruses = []
    for path in glob.glob(str(COHORT / "viruses" / "*.tsv")):
        for row in table(path):
            if not row.get("length"):
                continue
            quality = checkv[row["seq_name"]]
            ranks = row["taxonomy"].split(";")
            viruses.append({
                "name": row["seq_name"], "strategy": strategy[library(row["seq_name"])],
                "length": int(row["length"]), "hallmarks": int(row["n_hallmarks"]),
                "tier": quality["checkv_quality"], "method": quality["completeness_method"].split(" (")[0],
                "completeness": float(quality["completeness"]) if quality["completeness"] not in ("", "NA") else 0.0,
                "realm": ranks[1] if len(ranks) > 1 else "Unclassified",
                "phylum": ranks[3] if len(ranks) > 3 else "Unclassified",
            })
    good = lambda v: v["tier"] in ("Complete", "High-quality", "Medium-quality")
    rows = []
    for s in ("fixed_effort_trimmed", "reference_depleted"):
        for tier in ("Complete", "High-quality", "Medium-quality"):
            group = [v for v in viruses if v["strategy"] == s and v["tier"] == tier]
            lengths = sorted(v["length"] for v in group)
            rows.append([s, tier, len(group), lengths[len(lengths) // 2] if lengths else "",
                         sum(v["method"] == "DTR" for v in group), sum(v["hallmarks"] >= 1 for v in group),
                         sum(v["realm"] == "Unclassified" for v in group)])
    write("virus_checkv_tiers.tsv", ["strategy", "checkv_tier", "contigs", "median_length", "dtr_completeness",
                                     "with_hallmark", "unclassified"], rows)
    rules = {
        "genomad_default": lambda v: True,
        "checkv_medium_or_better": good,
        "hallmark": lambda v: v["hallmarks"] >= 1,
        "hallmark_and_5kb_or_50pct": lambda v: v["hallmarks"] >= 1 and (v["length"] >= 5000 or v["completeness"] >= 50),
        "hallmark_and_10kb_or_checkv_medium": lambda v: v["hallmarks"] >= 1 and (v["length"] >= 10000 or good(v)),
    }
    rows = []
    for label, rule in rules.items():
        kept = [v for v in viruses if rule(v)]
        rows.append([label, len(kept), len({library(v["name"]) for v in kept}),
                     sum(v["strategy"] == "fixed_effort_trimmed" for v in kept),
                     sum(v["strategy"] == "reference_depleted" for v in kept),
                     sum(v["realm"] == "Unclassified" for v in kept)])
    write("virus_rule_sensitivity.tsv", ["rule", "contigs", "libraries", "reference_free", "reference_depleted",
                                         "unclassified"], rows)
    selected = [v for v in viruses if rules["hallmark_and_5kb_or_50pct"](v)]
    write("virus_selected_phyla.tsv", ["phylum", "contigs"], Counter(v["phylum"] for v in selected).most_common())

    # Host-genome test: viral contigs aligned to the N. septata reference.
    by_name = {v["name"]: v for v in viruses}
    subsets = {"all": lambda v: True, "checkv_medium_or_better": good,
               "checkv_complete": lambda v: v["tier"] == "Complete", "hallmark": lambda v: v["hallmarks"] >= 1}
    rows = []
    for key in QUERY_SETS:
        members = [line[1:].split()[0] for line in open(OUT / f"{key}.fna") if line.startswith(">")]
        fractions = aligned_fraction(OUT / f"{key}.vs_septata.paf", {n: by_name[n]["length"] for n in members})
        for label, subset in subsets.items():
            group = [n for n in members if subset(by_name[n])]
            if group:
                rows.append([key, label, len(group), sum(fractions.get(n, 0) >= 0.5 for n in group)])
    write("virus_host_alignment.tsv", ["query_set", "subset", "contigs", "aligned_ge_50pct_to_N_septata"], rows)

    # The recurrent 5,000-bp Physalia DTR element.
    element = [v for v in viruses if v["length"] == 5000 and v["tier"] == "Complete"]
    write("physalia_5kb_element.tsv", ["contigs", "libraries", "studies", "physalis_reference_hits"],
          [[len(element), len({library(v["name"]) for v in element}),
            ",".join(sorted({v["name"].split("__")[0] for v in element})),
            sum(1 for _ in open(OUT / "physalia_5kb.vs_physalis.paf"))]])


parser = argparse.ArgumentParser()
parser.add_argument("step", choices=("prepare", "summarize"))
args = parser.parse_args()
prepare() if args.step == "prepare" else summarize()
