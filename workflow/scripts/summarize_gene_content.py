#!/usr/bin/env python3
"""Tabulate genome statistics, KEGG orthologs, and module completeness for the gene-content comparison.

KOs are assigned in two tiers: strict (KofamScan adaptive threshold) and relaxed (config kofam.relaxed);
relaxed tables are primary and strict ones are written alongside with a _strict suffix.
Module completeness comes from kegg-pathways-completeness (give_completeness) per genome.
"""
import argparse
import csv
import json
import os
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--table", required=True)
parser.add_argument("--checkm2", required=True)
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--genes", required=True, type=Path)
parser.add_argument("--kofam", required=True)
parser.add_argument("--ko-list", required=True)
parser.add_argument("--scratch", required=True, type=Path)
parser.add_argument("--genome-stats", required=True, type=Path)
parser.add_argument("--ko-matrix", required=True, type=Path)
parser.add_argument("--modules", required=True, type=Path)
parser.add_argument("--lineages", required=True, type=Path)
parser.add_argument("--focal", required=True, type=Path)
args = parser.parse_args()
settings = json.loads(Path(args.config).read_text())


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records, fields):
    temporary = path.with_name(path.name + ".tmp")
    path.parent.mkdir(parents=True, exist_ok=True)
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)
    os.replace(temporary, path)


lineage_of = {g: name for name, members in settings["lineages"].items() for g in members}
quality = {r["Name"]: r for r in rows(args.checkm2)}
definitions = {r["knum"]: r["definition"] for r in rows(args.ko_list)}
near = settings["near_complete"]

genomes = []
for record in rows(args.table):
    name = record["genome"]
    sequence = "".join(line.strip() for line in open(args.genomes / f"{name}.fa") if not line.startswith(">"))
    contigs = sum(1 for line in open(args.genomes / f"{name}.fa") if line.startswith(">"))
    genes = [line.split("\t") for line in open(args.genes / f"{name}.gff") if not line.startswith("#") and "\tCDS\t" in line]
    completeness = float(quality[name]["Completeness"])
    contamination = float(quality[name]["Contamination"])
    genomes.append({
        **record, "lineage": lineage_of.get(name, "this study" if name.startswith("MAGSP") else record["group"]),
        "completeness": completeness, "contamination": contamination,
        "included": completeness >= settings["min_completeness_included"]
        and contamination < settings["max_contamination_included"],
        "near_complete": completeness >= near["min_completeness"] and contamination < near["max_contamination"],
        "length": len(sequence), "contigs": contigs,
        "gc": round(sum(sequence.upper().count(b) for b in "GC") / len(sequence), 4),
        "genes": len(genes),
        "coding_density": round(sum(int(g[4]) - int(g[3]) + 1 for g in genes) / len(sequence), 4),
    })
by_name = {g["genome"]: g for g in genomes}

# KO assignment, one KO per protein. Strict: best threshold-passing hit. Relaxed (primary): a protein
# with no threshold-passing hit takes its best hit meeting the configured E-value and score fraction.
relaxed = settings["kofam"]["relaxed"]
hits = defaultdict(list)
for line in open(args.kofam):
    if line.startswith("#"):
        continue
    fields = line.rstrip("\n").split("\t")
    threshold = float(fields[3]) if fields[3] else None
    hits[fields[1]].append((fields[0] == "*", fields[2], float(fields[4]), threshold, float(fields[5])))
assignment = {"strict": {}, "relaxed": {}}
for protein, candidates in hits.items():
    passing = [h for h in candidates if h[0]]
    if passing:
        ko = max(passing, key=lambda h: h[2])[1]
        assignment["strict"][protein] = assignment["relaxed"][protein] = ko
        continue
    near_miss = [h for h in candidates if h[3] and h[4] <= relaxed["evalue_max"]
                 and h[2] >= relaxed["score_fraction_min"] * h[3]]
    if near_miss:
        assignment["relaxed"][protein] = max(near_miss, key=lambda h: h[2])[1]
counts = {}
for tier, assigned in assignment.items():
    counts[tier] = defaultdict(lambda: defaultdict(int))
    for protein, ko in assigned.items():
        counts[tier][protein.split("|", 1)[0]][ko] += 1
for g in genomes:
    for tier in counts:
        g[f"kos_{tier}"] = len(counts[tier][g["genome"]])
        g[f"proteins_with_ko_{tier}"] = sum(counts[tier][g["genome"]].values())

write(args.genome_stats, genomes, ["genome", "group", "lineage", "host_or_source", "classification", "completeness",
                                   "contamination", "included", "near_complete", "length", "contigs", "gc", "genes",
                                   "coding_density", "kos_relaxed", "proteins_with_ko_relaxed", "kos_strict",
                                   "proteins_with_ko_strict"])
included = [g["genome"] for g in genomes if g["included"]]
full = settings["module_complete_min"]
reference_sets = {
    "gtdb_family": [g for g in included if by_name[g]["near_complete"] and g.startswith("GTDB_")],
    "host_associated_external": [g for g in included if by_name[g]["near_complete"] and g.startswith("EXT_")],
}


def strict_path(path):
    return path.with_name(path.stem + "_strict" + path.suffix)


def tables(tier, ko_path, module_path, lineage_path):
    """KO matrix, per-genome module completeness, and the lineage view for one assignment tier."""
    tier_counts = counts[tier]
    all_kos = sorted({ko for g in included for ko in tier_counts[g]})
    write(ko_path, [{"ko": ko, "definition": definitions.get(ko, ""), **{g: tier_counts[g][ko] for g in included}}
                    for ko in all_kos], ["ko", "definition"] + included)

    args.scratch.mkdir(parents=True, exist_ok=True)
    modules, module_info = [], {}
    with tempfile.TemporaryDirectory(dir=args.scratch) as tmp:
        for name in included:
            listing = Path(tmp) / f"{name}.txt"
            listing.write_text(",".join(sorted(tier_counts[name])) + "\n")
            subprocess.run(["give_completeness", "-l", str(listing), "-o", str(Path(tmp) / name), "-r", name],
                           check=True, capture_output=True)
            for r in rows(Path(tmp) / name / f"{name}_pathways.tsv"):
                module_info[r["module_accession"]] = (r["pathway_name"], r["pathway_class"])
                modules.append({"genome": name, "module": r["module_accession"],
                                "completeness": float(r["completeness"]), "name": r["pathway_name"],
                                "class": r["pathway_class"], "missing_ko": r["missing_ko"]})
    write(module_path, modules, ["genome", "module", "completeness", "name", "class", "missing_ko"])

    # Lineage view: each siphonophore lineage against near-complete relatives.
    completeness = defaultdict(float)
    for m in modules:
        completeness[(m["genome"], m["module"])] = m["completeness"]
    records = []
    for module, (name, klass) in sorted(module_info.items()):
        record = {"module": module, "name": name, "class": klass}
        for lineage, members in settings["lineages"].items():
            members = [g for g in members if g in included]
            record[f"{lineage}: best"] = max((completeness[(g, module)] for g in members), default="")
            nc = [g for g in members if by_name[g]["near_complete"]]
            record[f"{lineage}: near_complete"] = ",".join(f"{g}={completeness[(g, module)]:g}" for g in nc)
        for label, members in reference_sets.items():
            record[f"{label}: fraction_complete"] = round(
                sum(completeness[(g, module)] >= full for g in members) / len(members), 3) if members else ""
            record[f"{label}: genomes"] = len(members)
        # Plan rule 3: complete in the lineage but rare in the family, or absent from all
        # near-complete lineage members but common in the family.
        family = record["gtdb_family: fraction_complete"]
        calls = []
        for lineage, members in settings["lineages"].items():
            members = [g for g in members if g in included]
            nc = [g for g in members if by_name[g]["near_complete"]]
            if family == "":
                continue
            if any(completeness[(g, module)] >= full for g in members) and family < 0.25:
                calls.append(f"{lineage}: complete, rare in family")
            if nc and all(completeness[(g, module)] == 0 for g in nc) and family >= 0.75:
                calls.append(f"{lineage}: not detected ({len(nc)} near-complete), common in family")
        record["distinctive"] = "; ".join(calls)
        records.append(record)
    fields = ["module", "name", "class"] + \
             [f"{l}: {k}" for l in settings["lineages"] for k in ("best", "near_complete")] + \
             [f"{l}: {k}" for l in reference_sets for k in ("fraction_complete", "genomes")] + ["distinctive"]
    write(lineage_path, records, fields)
    return len(module_info)


module_count = tables("relaxed", args.ko_matrix, args.modules, args.lineages)
tables("strict", strict_path(args.ko_matrix), strict_path(args.modules), strict_path(args.lineages))

# Focal functions: copies of each pre-specified KO per included genome under both tiers, with KOfam
# definitions as a check; the GTDB column is the fraction of near-complete family genomes carrying the KO.
focal = []
for function, kos in settings["focal_kos"].items():
    for ko in kos:
        for tier in ("relaxed", "strict"):
            family = reference_sets["gtdb_family"]
            focal.append({"function": function, "ko": ko, "definition": definitions.get(ko, "NOT IN KO LIST"),
                          "assignment": tier,
                          "gtdb_family_fraction": round(sum(counts[tier][g][ko] > 0 for g in family) / len(family), 3),
                          **{g: counts[tier][g][ko] for g in included}})
write(args.focal, focal, ["function", "ko", "definition", "assignment", "gtdb_family_fraction"] + included)
print(f"genomes={len(genomes)} included={len(included)} near_complete={sum(g['near_complete'] for g in genomes)} "
      f"proteins_with_ko strict={len(assignment['strict'])} relaxed={len(assignment['relaxed'])} modules={module_count}")
