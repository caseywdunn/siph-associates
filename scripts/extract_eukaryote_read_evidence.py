#!/usr/bin/env python3
"""Count filtered eukaryotic SSU read pairs per lineage from accepted phyloFlash read mappings.

Streams each library's nested phyloFlash bbmap.sam (reads mapped to SILVA 138.1)
and counts read pairs whose two mates both align over >= 100 bp to eukaryotic
references. Each pair gets the lineage its mates share (species labels dropped);
roll-up to reporting ranks happens in analysis. Pairs are binned by the lower of the
two mates' identity (>= 90, >= 95, >= 97%). Output is one long table,
data/results/eukaryote_gate/read_evidence.tsv. The archives are unchanged.
"""
from __future__ import annotations

import csv
import io
import re
import sys
import tarfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "data" / "results" / "screens" / "phyloflash"
OUT = ROOT / "data" / "results" / "eukaryote_gate" / "read_evidence.tsv"
MIN_ALIGNED = 100
BANDS = (97, 95, 90)
CIGAR = re.compile(r"(\d+)([=XIDMSHN])")


def alignment(fields):
    """Return (lineage, aligned_bases, identity_percent) for a mapped read, else None."""
    if fields[2] == "*" or int(fields[1]) & 0x4:
        return None
    reference = fields[2].split(" ", 1)
    lineage = reference[1] if len(reference) > 1 else ""
    matches = mismatches = aligned = 0
    for count, op in CIGAR.findall(fields[5]):
        n = int(count)
        if op == "=":
            matches += n; aligned += n
        elif op == "X":
            mismatches += n; aligned += n
        elif op in "ID":
            mismatches += n; aligned += n if op == "I" else 0
    return lineage, aligned, 100 * matches / (matches + mismatches) if matches + mismatches else 0


rows = Counter()
for archive in sorted(SOURCE.glob("*.tar.gz")):
    sample = archive.name.removesuffix(".tar.gz")
    with tarfile.open(archive) as outer:
        inner = next(m for m in outer.getmembers() if m.name.endswith(".phyloFlash.tar.gz"))
        with tarfile.open(fileobj=io.BytesIO(outer.extractfile(inner).read())) as nested:
            member = next((m for m in nested.getmembers() if m.name.endswith(".bbmap.sam")), None)
            if member is None:
                print(sample, "no bbmap.sam", file=sys.stderr)
                continue
            pending = {}
            for raw in io.TextIOWrapper(nested.extractfile(member), errors="replace"):
                if raw.startswith("@"):
                    continue
                fields = raw.rstrip("\n").split("\t")
                name = fields[0].split(" ")[0]
                hit = alignment(fields)
                if name not in pending:
                    pending[name] = hit
                    continue
                first = pending.pop(name)
                if not first or not hit or not first[0].startswith("Eukaryota") or not hit[0].startswith("Eukaryota"):
                    continue
                # A pair is assigned the lineage its two mates share, so mates hitting sister
                # genera (common for taxa not exactly represented in SILVA) are kept at the shared rank.
                shared = []
                for a, b in zip(first[0].split(";")[:-1], hit[0].split(";")[:-1]):
                    if a != b:
                        break
                    shared.append(a)
                lineage = ";".join(shared)
                if min(first[1], hit[1]) < MIN_ALIGNED:
                    continue
                identity = min(first[2], hit[2])
                band = next((b for b in BANDS if identity >= b), None)
                if band:
                    rows[(sample, lineage, band)] += 1
    print(sample, "done", file=sys.stderr)

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample_id", "lineage", "min_identity_band", "read_pairs"])
    for (sample, lineage, band), n in sorted(rows.items()):
        writer.writerow([sample, lineage, band, n])
print(f"rows={len(rows)}")
