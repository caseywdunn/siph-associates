#!/usr/bin/env python3
"""Select newly nominated genomes, then merge verified old and new reference QC."""

import argparse
from pathlib import Path

from assembly_completion import sha256, table, unique, write_table

ROOT = Path("data/results/post_expansion/2026-10-20/phase4_catalog/references")
OLD = Path("data/results/phase4_catalog/references")
NEW = ROOT / "new"


def select() -> None:
    nominated = table(ROOT / "nominated.tsv")
    previous = unique(table(OLD / "staged.tsv"), "genome_id")
    additions = [r for r in nominated if r["genome_id"] not in previous]
    write_table(NEW / "nominated.tsv", additions, list(nominated[0]))
    # Keep the established negative-control pool fixed across catalog versions.
    decoys = table(OLD / "decoy_candidates.tsv")
    write_table(ROOT / "decoy_candidates.tsv", decoys, list(decoys[0]))
    write_table(NEW / "decoy_candidates.tsv", [], list(decoys[0]))
    print(f"New reference genomes requiring QC: {len(additions)}")


def link(source: Path, destination: Path, expected: str) -> None:
    if sha256(source) != expected:
        raise ValueError(
            f"Reference sequence differs from its recorded identity: {source}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() == source.resolve():
            return
        destination.unlink()
    if destination.exists():
        raise FileExistsError(destination)
    destination.symlink_to(source.resolve())


def merge() -> None:
    merged = {}
    for name, key in [
        ("staged.tsv", "genome_id"),
        ("ncbi_status.tsv", "genome_id"),
        ("substitutes.tsv", "genome_id"),
        ("checkm2_quality.tsv", "Name"),
        ("oceandna_gtdbtk.bac120.summary.tsv", "user_genome"),
        ("oceandna_gtdbtk.ar53.summary.tsv", "user_genome"),
    ]:
        prior, new = table(OLD / name), table(NEW / name)
        combined = prior + new
        unique(combined, key)
        # Empty GTDB domain outputs may have only a user_genome header.
        fields = list((prior or new or [{key: ""}])[0])
        write_table(ROOT / name, combined, fields)
        merged[name] = combined
    for origin in (OLD, NEW):
        for row in table(origin / "staged.tsv"):
            link(
                origin / "genomes" / f"{row['genome_id']}.fa",
                ROOT / "genomes" / f"{row['genome_id']}.fa",
                row["fasta_sha256"],
            )
        for row in table(origin / "substitutes.tsv"):
            link(
                origin / "substitutes" / f"{row['substitute_accession']}.fa",
                ROOT / "substitutes" / f"{row['substitute_accession']}.fa",
                row["fasta_sha256"],
            )
    (ROOT / "substitutes").mkdir(exist_ok=True)
    required = {
        r["genome_id"]
        for r in table(ROOT / "nominated.tsv") + table(ROOT / "decoy_candidates.tsv")
    }
    assert required <= {r["genome_id"] for r in merged["staged.tsv"]}
    write_table(
        ROOT / "qc_sources.tsv",
        [
            {"source": str(p), "sha256": sha256(p)}
            for origin in (OLD, NEW)
            for name in merged
            for p in [origin / name]
        ],
        ["source", "sha256"],
    )
    print(f"Combined {len(merged['staged.tsv'])} staged reference/decoy genomes")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("select", "merge"))
    args = parser.parse_args()
    (select if args.stage == "select" else merge)()
