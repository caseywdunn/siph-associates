#!/usr/bin/env python3
"""Grade eukaryotic detections per library and reporting unit under config/eukaryote_gate.json.

Combines competitively remapped SSU read pairs (LCA-assigned), assembled SSUs from
phyloFlash and from the library's Phase-3 assembly (both LCA-classified against
the same reference), and NCBI verification of the assembled sequences.
"""
import argparse
import csv
import json
import os
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eukaryote_units import unit  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--config", required=True)
parser.add_argument("--reads", required=True)
parser.add_argument("--assembled", required=True, help="collected non-host assembled sequence table")
parser.add_argument("--verification", required=True)
parser.add_argument("--output", required=True, type=Path)
args = parser.parse_args()

gate = json.loads(Path(args.config).read_text())
presence, roles = gate["presence"], gate["role_field"]
band = presence["validated_reads"]["min_identity_band"]


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def named(lineage, identity):
    """Lineage reported to the depth its identity supports; SILVA species labels are never reported."""
    ranks = [r for r in lineage.split(";") if r]
    if identity >= 97:
        return ";".join(ranks), "genus"
    if identity >= 90:
        return ";".join(ranks[:-1]), "family_or_order"
    kind_unit = unit(lineage)
    return (kind_unit[1] if kind_unit else ";".join(ranks[:5])) + " (novel lineage)", "class_or_phylum"


def role(lineage, ncbi_hits):
    """Verified NCBI overrides first (e.g. the Paramoeba endosymbiont), then the lineage lists."""
    for override in roles.get("verified_overrides", []):
        if any(override["ncbi_top_hit_contains"] in h for h in ncbi_hits):
            return override["role"]
    ranks = set(lineage.split(";"))
    for label in ("parasite", "prey"):
        if ranks & set(roles[label]):
            return label
    return "unassigned"


evidence = defaultdict(lambda: {"read_pairs": 0, "read_lineages": Counter(), "assembled": []})
kinds = {}
for r in rows(args.reads):
    u = unit(r["lineage"])
    if u and int(r["min_identity_band"]) >= band:
        kinds[u[1]] = u[0]
        if u[0] in ("metazoan", "non_metazoan", "human"):
            entry = evidence[(r["sample_id"], u[1])]
            entry["read_pairs"] += int(r["read_pairs"])
            entry["read_lineages"][r["lineage"]] += int(r["read_pairs"])
verified = defaultdict(list)
for r in rows(args.verification):
    if r["rank"] == "1":
        verified[r["sequence_id"]].append(r)
for r in rows(args.assembled):
    u = unit(r["lineage"])
    if u:
        kinds[u[1]] = u[0]
        evidence[(r["sample_id"], u[1])]["assembled"].append(r)

records = []
for (sample, reporting_unit), e in sorted(evidence.items()):
    reads_ok = e["read_pairs"] >= presence["validated_reads"]["min_read_pairs"]
    trace = e["read_pairs"] >= presence["trace_reads"]["min_read_pairs"]
    assembled = e["assembled"]
    if reads_ok and assembled:
        grade = "high_confidence"
    elif reads_ok or assembled:
        grade = "validated"
    elif trace:
        grade = "trace"
    else:
        continue
    if assembled:
        best = max(assembled, key=lambda a: float(a["identity"]))
        lineage, depth = named(best["lineage"], float(best["identity"]))
        identity = float(best["identity"])
    else:
        top = e["read_lineages"].most_common(1)[0][0]
        lineage, depth, identity = named(top, band)[0], named(top, band)[1], band
    top_hits = [verified[a["sequence_id"]][0] for a in assembled if verified.get(a["sequence_id"])]
    records.append({
        "sample_id": sample, "reporting_unit": reporting_unit, "kind": kinds[reporting_unit], "grade": grade,
        "read_pairs_ge97": e["read_pairs"], "assembled_sequences": len(assembled),
        "assembled_sources": ",".join(sorted({a["source"] for a in assembled})),
        "best_identity": round(identity, 2), "named_lineage": lineage, "naming_depth": depth,
        "role": role(lineage, [h["organism"] + " " + h.get("title", "") for h in top_hits]), "contamination": str(kinds[reporting_unit] == "human").lower(),
        "ncbi_top_hit": "; ".join(f"{h['organism']} ({h['subject']}, {h['identity']}%)" for h in top_hits[:2]),
    })

fields = list(records[0]) if records else ["sample_id"]
temporary = args.output.with_name(args.output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
os.replace(temporary, args.output)
print("grades:", dict(Counter(r["grade"] for r in records)),
      "roles:", dict(Counter(r["role"] for r in records if r["grade"] != "trace")))
