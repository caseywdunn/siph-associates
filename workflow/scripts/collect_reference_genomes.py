#!/usr/bin/env python3
"""Stage nominated and decoy-candidate genomes as uncompressed FASTA with checksums."""
import argparse
import csv
import gzip
import hashlib
import os
import shutil
import tarfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--references", required=True)
parser.add_argument("--decoys", required=True)
parser.add_argument("--gtdb-genomes", required=True, type=Path)
parser.add_argument("--oceandna-archive", required=True)
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--table", required=True, type=Path)
parser.add_argument("--oceandna-batch", required=True, type=Path)
args = parser.parse_args()


def rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def gtdb_path(accession):
    # GCF_000354175.2 -> GCF/000/354/175/GCF_000354175.2_genomic.fna.gz
    prefix, digits = accession[:3], accession[4:13]
    return args.gtdb_genomes / prefix / digits[0:3] / digits[3:6] / digits[6:9] / f"{accession}_genomic.fna.gz"


args.genomes.mkdir(parents=True)
staged = []
wanted = {}
for role, table in (("nominated", args.references), ("decoy_candidate", args.decoys)):
    for row in rows(table):
        source = "OceanDNA" if row["genome_id"].startswith("OceanDNA") else "GTDB"
        wanted[row["genome_id"]] = (role, source)

for genome_id, (role, source) in wanted.items():
    if source != "GTDB":
        continue
    origin = gtdb_path(genome_id)
    target = args.genomes / f"{genome_id}.fa"
    with gzip.open(origin, "rb") as src, open(target, "wb") as dst:
        shutil.copyfileobj(src, dst)
    staged.append((genome_id, role, source, str(origin), sha256(origin), target))

# One sequential pass over the OceanDNA archive extracts the nominated members.
oceandna = {g for g, (_, s) in wanted.items() if s == "OceanDNA"}
with tarfile.open(args.oceandna_archive) as archive:
    for member in archive:
        name = Path(member.name).name
        genome_id = name.removesuffix(".fa.gz")
        if member.isfile() and genome_id in oceandna:
            target = args.genomes / f"{genome_id}.fa"
            with gzip.open(archive.extractfile(member), "rb") as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst)
            staged.append((genome_id, wanted[genome_id][0], "OceanDNA", f"{args.oceandna_archive}:{member.name}",
                           "", target))
            oceandna.discard(genome_id)
if oceandna:
    raise SystemExit(f"OceanDNA genomes missing from the archive: {sorted(oceandna)}")

temporary = args.table.with_name(args.table.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["genome_id", "role", "source", "origin", "origin_sha256", "fasta_sha256"])
    for genome_id, role, source, origin, origin_sha, target in staged:
        writer.writerow([genome_id, role, source, origin, origin_sha, sha256(target)])
os.replace(temporary, args.table)
with args.oceandna_batch.open("w") as handle:
    for genome_id, _, source, _, _, target in staged:
        if source == "OceanDNA":
            handle.write(f"{target}\t{genome_id}\n")
print(f"staged {len(staged)} genomes")
