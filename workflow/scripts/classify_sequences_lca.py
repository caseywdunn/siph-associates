#!/usr/bin/env python3
"""Assign each query sequence to the LCA of its near-tied hits (single sequences, not pairs).

Reads minimap2 SAM on stdin (secondaries kept). For each query, hits aligned
over >= --min-aligned bp contribute their lineages (species label dropped);
identity is the best hit's. Used for assembled 18S genes.
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--taxonomy", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--min-aligned", type=int, default=100)
args = parser.parse_args()
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
NM = re.compile(r"\tNM:i:(\d+)")

lineages = {}
for line in open(args.taxonomy):
    ref, lineage = line.rstrip("\n").split("\t", 1)
    ranks = lineage.split(";")
    lineages[ref] = tuple(ranks[:-1] if len(ranks) > 1 else ranks)

hits = defaultdict(lambda: [set(), 0.0, 0])
for line in sys.stdin:
    if line.startswith("@"):
        continue
    f = line.split("\t", 11)
    if int(f[1]) & 0x4 or f[2] == "*":
        continue
    aligned = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in "M=XI")
    if aligned < args.min_aligned:
        continue
    match = NM.search(line)
    identity = 100 * (aligned - (int(match.group(1)) if match else 0)) / aligned
    record = hits[f[0]]
    record[0].add(lineages.get(f[2], ("unknown",)))
    if identity > record[1]:
        record[1], record[2] = identity, aligned

with open(args.output + ".tmp", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample_id", "sequence_id", "lineage", "identity", "aligned_bases"])
    for query, (groups, identity, aligned) in sorted(hits.items()):
        shared = []
        for ranks in zip(*groups):
            if len(set(ranks)) != 1:
                break
            shared.append(ranks[0])
        writer.writerow([args.sample, query, ";".join(shared), round(identity, 2), aligned])
os.replace(args.output + ".tmp", args.output)
print(f"{args.sample}: classified {len(hits)} sequences", file=sys.stderr)
