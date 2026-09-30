#!/usr/bin/env python3
"""Assign SSU read pairs to the lowest common ancestor of their near-tied hits.

Reads SAM from minimap2 (-N 200 -p 0.99, secondaries kept) on stdin. For each
mate, hits aligned over >= 100 bp are kept; the mate's lineage set is the
lineages of those hits (species labels dropped). A pair's lineage is the LCA of
both mates' sets, so reads from conserved regions resolve only to broad ranks
and host reads resolve to the host references. Pairs are binned by the lower
of the two mates' best identity (>= 90, 95, 97, 99%).
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--taxonomy", required=True, help="reference id <tab> lineage")
parser.add_argument("--sample", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--min-aligned", type=int, default=100)
args = parser.parse_args()
BANDS = (99, 97, 95, 90)
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
NM = re.compile(r"\tNM:i:(\d+)")

lineages = {}
with open(args.taxonomy) as handle:
    for line in handle:
        ref, lineage = line.rstrip("\n").split("\t", 1)
        ranks = lineage.split(";")
        lineages[ref] = tuple(ranks[:-1] if len(ranks) > 1 else ranks)

mates = defaultdict(lambda: [set(), set(), 0.0, 0.0])  # lineages mate1, mate2, best identity mate1, mate2
for line in sys.stdin:
    if line.startswith("@"):
        continue
    f = line.split("\t", 12)
    flag = int(f[1])
    if flag & 0x4 or f[2] == "*":
        continue
    aligned = sum(int(n) for n, op in CIGAR.findall(f[5]) if op in "M=XI")
    if aligned < args.min_aligned:
        continue
    match = NM.search(line)
    nm = int(match.group(1)) if match else 0
    identity = 100 * (aligned - nm) / aligned
    slot = 0 if flag & 0x40 else 1
    record = mates[f[0]]
    record[slot].add(lineages.get(f[2], ("unknown",)))
    record[2 + slot] = max(record[2 + slot], identity)


def lca(sets):
    groups = [g for s in sets for g in s]
    shared = []
    for ranks in zip(*groups):
        if len(set(ranks)) != 1:
            break
        shared.append(ranks[0])
    return ";".join(shared)


counts = Counter()
for first, second, id1, id2 in mates.values():
    if not first or not second:
        continue
    identity = min(id1, id2)
    band = next((b for b in BANDS if identity >= b), None)
    if band:
        counts[(lca([first, second]), band)] += 1

with open(args.output + ".tmp", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample_id", "lineage", "min_identity_band", "read_pairs"])
    for (lineage, band), n in sorted(counts.items()):
        writer.writerow([args.sample, lineage, band, n])
os.replace(args.output + ".tmp", args.output)
print(f"{args.sample}: pairs classified={sum(counts.values())}", file=sys.stderr)
