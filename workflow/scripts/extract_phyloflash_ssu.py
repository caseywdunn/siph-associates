#!/usr/bin/env python3
"""Extract one library's SSU reads and assembled SSUs from its accepted phyloFlash archive.

The archive is verified against the SHA-256 recorded in its Phase-2 provenance
before anything is read; it is passed as a path rather than a workflow input so
the known trim-provenance mtime cascade cannot reopen the Phase-2 screens.
"""
import argparse
import gzip
import hashlib
import io
import json
import os
import shutil
import tarfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--sample", required=True)
parser.add_argument("--provenance", required=True)
parser.add_argument("--r1", required=True, type=Path)
parser.add_argument("--r2", required=True, type=Path)
parser.add_argument("--assembled", required=True, type=Path)
parser.add_argument("--classification", required=True, type=Path)
args = parser.parse_args()

record = json.loads(Path(args.provenance).read_text())["archive"]
digest = hashlib.sha256()
with open(record["path"], "rb") as handle:
    for block in iter(lambda: handle.read(1 << 20), b""):
        digest.update(block)
if digest.hexdigest() != record["sha256"]:
    raise SystemExit(f"archive checksum differs from Phase-2 provenance: {record['path']}")

targets = {".SSU.1.fq": args.r1, ".SSU.2.fq": args.r2, ".all.final.fasta": args.assembled,
           ".phyloFlash.extractedSSUclassifications.csv": args.classification}
found = set()
with tarfile.open(record["path"]) as outer:
    inner = next(m for m in outer.getmembers() if m.name.endswith(".phyloFlash.tar.gz"))
    with tarfile.open(fileobj=io.BytesIO(outer.extractfile(inner).read())) as nested:
        for member in nested.getmembers():
            name = Path(member.name).name
            for suffix, target in targets.items():
                if member.isfile() and name.endswith(suffix) and suffix not in found:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    temporary = target.with_name(target.name + ".tmp")
                    source = nested.extractfile(member)
                    if target.suffix == ".gz":
                        with gzip.open(temporary, "wb") as out:
                            shutil.copyfileobj(source, out)
                    else:
                        temporary.write_bytes(source.read())
                    os.replace(temporary, target)
                    found.add(suffix)
# phyloFlash writes no assembled SSU when SPAdes finds no rRNA contig; record that as empty.
for suffix in (".all.final.fasta", ".phyloFlash.extractedSSUclassifications.csv"):
    if suffix not in found:
        targets[suffix].parent.mkdir(parents=True, exist_ok=True)
        targets[suffix].write_text("" if suffix.endswith("fasta") else "OTU,read_cov,coverage,dbHit,taxonomy,%id,alnlen,evalue\n")
missing = [s for s in (".SSU.1.fq", ".SSU.2.fq") if s not in found]
if missing:
    raise SystemExit(f"SSU reads missing from archive: {missing}")
print(f"{args.sample}: extracted {sorted(found)}")
