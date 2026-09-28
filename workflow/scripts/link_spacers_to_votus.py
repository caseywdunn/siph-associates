#!/usr/bin/env python3
"""Keep full-length spacer matches with at most one mismatch and summarize host links."""
import argparse
import csv
import json
import os
from collections import defaultdict
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--blast", required=True, help="outfmt '6 std qlen slen'")
parser.add_argument("--manifest", required=True)
parser.add_argument("--config", required=True)
parser.add_argument("--matches", required=True, type=Path)
parser.add_argument("--links", required=True, type=Path)
args = parser.parse_args()

rule = json.loads(Path(args.config).read_text())["viruses"]["crispr_host_linkage"]
with open(args.manifest, newline="") as handle:
    taxonomy = {r["catalog_id"]: r["taxonomy"] for r in csv.DictReader(handle, delimiter="\t") if r["catalog_id"]}

kept = []
for line in open(args.blast):
    f = line.rstrip("\n").split("\t")
    spacer, target = f[0], f[1]
    length, mismatch, gaps, qstart, qend, qlen = int(f[3]), int(f[4]), int(f[5]), int(f[6]), int(f[7]), int(f[12])
    full = qstart == 1 and qend == qlen and length >= qlen
    if (full or not rule["require_full_spacer_length"]) and mismatch + gaps <= rule["max_mismatches"]:
        kept.append((spacer, target, mismatch + gaps, qlen))

links = defaultdict(set)
with args.matches.with_name(args.matches.name + ".tmp").open("w") as out:
    out.write("spacer_id\tvotu_id\tmismatches\tspacer_length\thost_catalog_id\n")
    for spacer, target, mismatches, qlen in kept:
        votu = target.split()[0]
        host = spacer.split("|")[0]
        out.write(f"{spacer}\t{votu}\t{mismatches}\t{qlen}\t{host}\n")
        links[(votu, host)].add(spacer)
os.replace(args.matches.with_name(args.matches.name + ".tmp"), args.matches)
with args.links.with_name(args.links.name + ".tmp").open("w") as out:
    out.write("votu_id\thost_catalog_id\tspacers\thost_taxonomy\n")
    for (votu, host), spacers in sorted(links.items()):
        out.write(f"{votu}\t{host}\t{len(spacers)}\t{taxonomy.get(host, '')}\n")
os.replace(args.links.with_name(args.links.name + ".tmp"), args.links)
print(f"spacer matches={len(kept)} host links={len(links)} vOTUs linked={len({v for v, _ in links})}")
