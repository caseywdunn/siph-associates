#!/usr/bin/env python3
"""Separate associate viruses from endogenous candidates using host alignment and taxonomy."""
import argparse
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--included", required=True)
parser.add_argument("--fasta", required=True)
parser.add_argument("--paf", required=True, nargs="+", help="ROUTE=PATH alignments to host references")
parser.add_argument("--samples", required=True)
parser.add_argument("--manifest", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--table", required=True, type=Path)
parser.add_argument("--associate", required=True, type=Path)
parser.add_argument("--endogenous", required=True, type=Path)
args = parser.parse_args()

rules = json.loads(Path(args.config).read_text())["viruses"]["endogenous_candidates"]
with open(args.samples, newline="") as handle:
    route = {row["sample_id"]: row["host_route"] for row in csv.DictReader(handle, delimiter="\t")}
with open(args.manifest, newline="") as handle:
    species = {f"{row['study']}__{row['library_id'].split(':')[-1]}": row["species_current"]
               for row in csv.DictReader(handle)}


def tested_reference(sample):
    """The host reference that can reveal host origin: the library's own, or the congener for Nanomia."""
    if route[sample] != "none":
        return route[sample]
    return "N_septata" if species[sample].startswith("Nanomia") else ""


def aligned_fraction(paf):
    spans = defaultdict(list)
    lengths = {}
    for line in open(paf):
        fields = line.split("\t")
        spans[fields[0]].append((int(fields[2]), int(fields[3])))
        lengths[fields[0]] = int(fields[1])
    fractions = {}
    for query, intervals in spans.items():
        intervals.sort()
        total, (start, end) = 0, intervals[0]
        for s, e in intervals[1:]:
            if s > end:
                total, start, end = total + end - start, s, e
            else:
                end = max(end, e)
        fractions[query] = (total + end - start) / lengths[query]
    return fractions


fractions = {}
for item in args.paf:
    reference, path = item.split("=", 1)
    fractions[reference] = aligned_fraction(path)

with open(args.included, newline="") as handle:
    viruses = list(csv.DictReader(handle, delimiter="\t"))
for virus in viruses:
    reference = tested_reference(virus["sample_id"])
    fraction = fractions[reference].get(virus["contig_id"], 0.0) if reference else None
    reasons = []
    if fraction is not None and fraction >= rules["min_host_aligned_fraction"]:
        reasons.append(f"host_aligned_{reference}")
    if virus["phylum"] in rules["flag_phyla"]:
        reasons.append(f"endogenous_prone_{virus['phylum']}")
    virus["host_reference_tested"] = reference or "none_available"
    virus["host_aligned_fraction"] = "" if fraction is None else f"{fraction:.3f}"
    virus["virus_class"] = "endogenous_candidate" if reasons else "associate"
    virus["endogenous_reasons"] = ",".join(reasons)

classes = {virus["contig_id"]: virus["virus_class"] for virus in viruses}
outputs = {"associate": args.associate, "endogenous_candidate": args.endogenous}
handles = {key: path.with_name(path.name + ".tmp").open("w") for key, path in outputs.items()}
target = None
for line in open(args.fasta):
    if line.startswith(">"):
        target = handles[classes[line[1:].split()[0]]]
    target.write(line)
for key, handle in handles.items():
    handle.close()
    os.replace(outputs[key].with_name(outputs[key].name + ".tmp"), outputs[key])

temporary = args.table.with_name(args.table.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(viruses[0]) if viruses else ["contig_id"],
                            delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(viruses)
os.replace(temporary, args.table)
print(f"associate={sum(v == 'associate' for v in classes.values())} "
      f"endogenous_candidate={sum(v == 'endogenous_candidate' for v in classes.values())}")
