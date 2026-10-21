#!/usr/bin/env python3
"""Prepare candidate contigs and retain only CheckV viral regions for review."""

import json
from collections import defaultdict
from pathlib import Path

from assembly_completion import sha256, table, unique, write_table
from Bio import SeqIO

ROOT = Path("data/results/post_expansion/2026-10-20")
OUT = ROOT / "qc/curation"


def viral_regions(row: dict) -> list[tuple[int, int]]:
    """Return CheckV's one-based inclusive viral intervals, preserving splits."""
    length = int(row["contig_length"])
    if row["provirus"] == "No":
        return [(1, length)]
    regions, last = [], 0
    for kind, coords in zip(
        row["region_types"].split(","), row["region_coords_bp"].split(","), strict=True
    ):
        start, end = map(int, coords.split("-"))
        if start != last + 1 or not start <= end <= length:
            raise ValueError(f"Invalid CheckV region partition: {row['contig_id']}")
        if kind == "viral":
            regions.append((start, end))
        elif kind != "host":
            raise ValueError(f"Unknown CheckV region: {kind}")
        last = end
    if last != length or not regions:
        raise ValueError(f"Incomplete CheckV region partition: {row['contig_id']}")
    return regions


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    bins = table("figures/catalog_qc/domain_review.tsv")
    selected = [r for r in bins if int(r["genes_Bacteria"]) > int(r["genes_Eukaryota"])]
    seen = set()
    with (OUT / "bacterial_candidates.fna").open("w") as out:
        for row in selected:
            if sha256(row["fasta"]) != row["sha256"]:
                raise ValueError(f"Candidate changed: {row['mag_id']}")
            for record in SeqIO.parse(row["fasta"], "fasta"):
                if record.id in seen:
                    raise ValueError(f"Repeated candidate contig: {record.id}")
                seen.add(record.id)
                SeqIO.write(record, out, "fasta")
    controls = json.loads(Path("config/mag_domain_audit.json").read_text())["bins"]
    with (OUT / "perkinsela_controls.fna").open("w") as out:
        for row in controls:
            if sha256(row["fasta"]) != row["sha256"]:
                raise ValueError(f"Control changed: {row['bin_id']}")
            SeqIO.write(SeqIO.parse(row["fasta"], "fasta"), out, "fasta")

    source = ROOT / "phase3_catalog/viruses"
    qc = ROOT / "qc/viruses"
    candidates = unique(table(source / "included.tsv"), "contig_id")
    checkv = unique(table(qc / "checkv/contamination.tsv"), "contig_id")
    genes = defaultdict(list)
    for row in table(qc / "genomad/included_annotate/included_genes.tsv"):
        if row["virus_hallmark"] == "1":
            genes[row["gene"].rsplit("_", 1)[0]].append(row)
    records = []
    with (OUT / "viral_regions.fna").open("w") as out:
        for record in SeqIO.parse(source / "included.fna", "fasta"):
            if len(record) != int(checkv[record.id]["contig_length"]):
                raise ValueError(f"Viral sequence length mismatch: {record.id}")
            for start, end in viral_regions(checkv[record.id]):
                trimmed = (start, end) != (1, len(record))
                region_id = f"{record.id}|qc_{start}_{end}" if trimmed else record.id
                hallmarks = [
                    g
                    for g in genes[record.id]
                    if start <= int(g["start"]) <= int(g["end"]) <= end
                ]
                row = dict(
                    candidates[record.id],
                    contig_id=region_id,
                    source_contig_id=record.id,
                    source_start=start,
                    source_end=end,
                    source_length=len(record),
                    length=end - start + 1,
                    genomad_hallmarks=len(hallmarks),
                    hallmark_functions="; ".join(
                        g["annotation_description"] for g in hallmarks
                    ),
                    cellular_flanks_removed=str(trimmed).lower(),
                )
                records.append(row)
                region = record[start - 1 : end]
                region.id, region.description = region_id, ""
                SeqIO.write(region, out, "fasta")
    write_table(OUT / "viral_regions.tsv", records, list(records[0]))
    print(
        json.dumps(
            {
                "bacterial_candidates": len(selected),
                "bacterial_contigs": len(seen),
                "viral_source_contigs": len(candidates),
                "viral_regions": len(records),
                "trimmed_regions": sum(
                    r["cellular_flanks_removed"] == "true" for r in records
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
