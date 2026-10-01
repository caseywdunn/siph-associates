#!/usr/bin/env python3
"""Download NCBI genome assemblies listed in a table (accession, ftp) and verify NCBI's MD5."""
import argparse
import csv
import gzip
import hashlib
import os
import shutil
import time
import urllib.request
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--table", required=True)
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--manifest", required=True, type=Path)
args = parser.parse_args()


def fetch(url, path):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300) as response, open(path, "wb") as out:
                shutil.copyfileobj(response, out)
            return
        except Exception:
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"download failed: {url}")


def digest(path, algorithm):
    h = hashlib.new(algorithm)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


args.genomes.mkdir(parents=True)
records = []
with open(args.table, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        base = row["ftp"].rstrip("/")
        name = f"{base.rsplit('/', 1)[1]}_genomic.fna.gz"
        archive = args.genomes / name
        fetch(f"{base}/{name}", archive)
        checksums = args.genomes / f"{row['accession']}.md5checksums.txt"
        fetch(f"{base}/md5checksums.txt", checksums)
        expected = next(line.split()[0] for line in checksums.read_text().splitlines() if line.endswith(f"/{name}"))
        if digest(archive, "md5") != expected:
            raise SystemExit(f"MD5 mismatch for {name}")
        fasta = args.genomes / f"{row['accession']}.fa"
        with gzip.open(archive, "rb") as src, open(fasta, "wb") as dst:
            shutil.copyfileobj(src, dst)
        archive.unlink()
        checksums.unlink()
        records.append([row["accession"], f"{base}/{name}", expected, digest(fasta, "sha256")])

temporary = args.manifest.with_name(args.manifest.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["accession", "url", "ncbi_md5", "fasta_sha256"])
    writer.writerows(records)
os.replace(temporary, args.manifest)
print(f"fetched {len(records)} genomes")
