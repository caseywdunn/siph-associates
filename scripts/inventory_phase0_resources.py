#!/usr/bin/env python3
"""Inventory and checksum the reference/database resources locked at Phase 0."""
from __future__ import annotations

import csv
import hashlib
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "metadata" / "input_resources.tsv"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


resources = [
    dict(resource_id="host_physalia", category="host_reference",
         version="GCA_041430235.2",
         path="/gpfs/ycga/work/dunn/cwd7/databases/refgenomes/physalia_physalis/GCA_041430235.2_Physalia_TX2017-38_primary_02_genomic.fna",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/refgenomes/physalia_physalis/GCA_041430235.2_Physalia_TX2017-38_primary_02_genomic.fna.gz",
         notes="complete NCBI genomic download; whole-genome subtraction target"),
    dict(resource_id="host_nanomia_septata", category="host_reference",
         version="GCA_048301705.1",
         path="/gpfs/ycga/work/dunn/cwd7/databases/refgenomes/nanomia_septata/GCA_048301705.1_ASM4830170v1_genomic.fna",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/refgenomes/nanomia_septata/GCA_048301705.1_ASM4830170v1_genomic.fna.gz",
         notes="complete NCBI genomic download; whole-genome subtraction target"),
    dict(resource_id="kraken2_standard", category="screen_database",
         version="k2_standard_20240904", path="/gpfs/gibbs/data/db/kraken2/standard",
         anchor="/gpfs/gibbs/data/db/kraken2/k2_standard_20240904/inspect.txt",
         notes="81 GB standard database; anchor is Kraken inspection report"),
    dict(resource_id="sylph_gtdb", category="screen_database", version="GTDB r220 c200",
         path="/gpfs/ycga/work/dunn/cwd7/databases/sylph_gtdb_r220/gtdb-r220-c200-dbv1.syldb",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/sylph_gtdb_r220/gtdb-r220-c200-dbv1.syldb",
         notes="full database file checksummed"),
    dict(resource_id="sylph_oceandna", category="screen_database", version="OceanDNA c200 v0.3",
         path="/gpfs/ycga/work/dunn/cwd7/databases/sylph_gtdb_r220/OceanDNA-c200-v0.3.syldb",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/sylph_gtdb_r220/OceanDNA-c200-v0.3.syldb",
         notes="full database file checksummed"),
    dict(resource_id="phyloflash_silva", category="screen_database", version="SILVA 138.1 NR99",
         path="/gpfs/ycga/work/dunn/cwd7/databases/phyloflash_silva/138.1",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/phyloflash_silva/138.1/SILVA_SSU.noLSU.masked.trimmed.NR99.fixed.fasta",
         notes="local phyloFlash build; source FASTA checksummed"),
    dict(resource_id="gtdbtk", category="taxonomy_database", version="r220",
         path="/gpfs/ycga/work/dunn/cwd7/databases/gtdbtk_r220/release220",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/gtdbtk_r220/release220/metadata/metadata.txt",
         notes="VERSION_DATA=r220 metadata anchor"),
    dict(resource_id="checkm2", category="quality_database", version="local pilot install",
         path="/gpfs/ycga/work/dunn/cwd7/databases/checkm2_db/CheckM2_database",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd",
         notes="full DIAMOND database checksummed"),
    dict(resource_id="genomad", category="virus_database", version="1.9",
         path="/gpfs/ycga/work/dunn/cwd7/databases/genomad_db/genomad_db",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/genomad_db/genomad_db/genomad_db",
         notes="full main MMseqs database component checksummed"),
    dict(resource_id="checkv", category="virus_database", version="1.5",
         path="/gpfs/ycga/work/dunn/cwd7/databases/checkv_db/checkv-db-v1.5",
         anchor="/gpfs/ycga/work/dunn/cwd7/databases/checkv_db/checkv-db-v1.5/genome_db/checkv_info.tsv",
         notes="versioned database metadata anchor"),
]

columns = ["resource_id", "category", "version", "path", "status",
           "path_size_bytes", "checksum_anchor", "anchor_size_bytes",
           "anchor_sha256", "notes"]
rows = []
for resource in resources:
    path, anchor = Path(resource["path"]), Path(resource["anchor"])
    status = "ready" if path.exists() and anchor.is_file() else "missing"
    rows.append(dict(
        resource_id=resource["resource_id"], category=resource["category"],
        version=resource["version"], path=path, status=status,
        path_size_bytes=(path.stat().st_size if path.is_file() else ""),
        checksum_anchor=anchor,
        anchor_size_bytes=(anchor.stat().st_size if anchor.is_file() else ""),
        anchor_sha256=(digest(anchor) if anchor.is_file() else ""),
        notes=resource["notes"],
    ))

OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
print(f"wrote {len(rows)} resources -> {OUT}")
if any(row["status"] != "ready" for row in rows):
    raise SystemExit("one or more Phase-0 resources are missing")
