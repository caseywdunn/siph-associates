#!/usr/bin/env python3
"""Extract eukaryotic SSU evidence from accepted Phase-2 phyloFlash archives.

For every library, copy the assembled full-length SSU sequences, their SILVA
best-hit classifications, and the read-level NTU table out of the nested
phyloFlash archives into data/results/eukaryote_gate/evidence/<sample>/.
The archives themselves are unchanged.
"""
from __future__ import annotations

import io
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "data" / "results" / "screens" / "phyloflash"
OUT = ROOT / "data" / "results" / "eukaryote_gate" / "evidence"
WANTED = (".phyloFlash.extractedSSUclassifications.csv", ".all.final.fasta", ".phyloFlash.NTUfull_abundance.csv")

for archive in sorted(SOURCE.glob("*.tar.gz")):
    sample = archive.name.removesuffix(".tar.gz")
    target = OUT / sample
    if target.is_dir() and all((target / f"{sample}{suffix}").is_file() for suffix in WANTED):
        continue
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as outer:
        inner_member = next(m for m in outer.getmembers() if m.name.endswith(".phyloFlash.tar.gz"))
        with tarfile.open(fileobj=io.BytesIO(outer.extractfile(inner_member).read())) as inner:
            for member in inner.getmembers():
                name = Path(member.name).name
                if member.isfile() and name.endswith(WANTED):
                    (target / name).write_bytes(inner.extractfile(member).read())
    missing = [s for s in WANTED if not (target / f"{sample}{s}").is_file()]
    print(sample, "missing:" if missing else "ok", *missing)
