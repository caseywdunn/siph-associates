#!/usr/bin/env python3
"""Select discordant genes and their contigs for broader reference adjudication."""

import csv
import json
from collections import defaultdict
from pathlib import Path

from assembly_completion import sha256, table, unique, write_table
from Bio import SeqIO


def main() -> None:
    config = json.loads(Path("config/catalog_nr_review.json").read_text())
    root, source = Path(config["run"]), Path(config["screen_root"])
    bins = unique(table("figures/catalog_qc/domain_review.tsv"), "mag_id")
    bacterial = {
        k
        for k, r in bins.items()
        if int(r["genes_Bacteria"]) > int(r["genes_Eukaryota"])
    }
    contigs = table("data/results/post_expansion/2026-10-20/qc/contig_evidence.tsv")
    flagged = {
        r["contig_id"]
        for r in contigs
        if r["mag_id"] in bacterial
        and int(r["genes_Eukaryota"]) >= 3
        and int(r["genes_Eukaryota"]) > int(r["genes_Bacteria"])
    }
    chosen, euk = {}, defaultdict(list)
    with (source / "review/genes_swissprot.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            reasons = []
            if row["bin_id"] in bacterial and row["assignment"] in (
                "Eukaryota",
                "Archaea",
            ):
                reasons.append("nonbacterial_gene_in_bacterial_candidate")
            if row["contig_id"] in flagged:
                reasons.append("all_genes_on_eukaryote_enriched_contig")
            if reasons:
                chosen[row["query_id"]] = ";".join(reasons)
            if row["bin_id"] not in bacterial and row["assignment"] == "Eukaryota":
                euk[row["bin_id"]].append(row)
    # Twenty longest supported eukaryotic proteins per excluded candidate clarify
    # taxonomic identity. They do not replace the whole-genome domain evidence.
    for rows in euk.values():
        for row in sorted(rows, key=lambda r: (-int(r["length"]), r["query_id"]))[:20]:
            chosen[row["query_id"]] = "eukaryotic_candidate_identity_control"
    genes = [r for r in table(source / "genes.tsv") if r["query_id"] in chosen]
    assert len(genes) == len(chosen)
    root.mkdir(parents=True, exist_ok=True)
    with (root / "genes.faa").open("w") as out:
        written = set()
        for r in SeqIO.parse(source / "genes.faa", "fasta"):
            if r.id in chosen:
                SeqIO.write(r, out, "fasta")
                written.add(r.id)
    assert written == set(chosen)
    write_table(root / "genes.tsv", genes, list(genes[0]))
    write_table(
        root / "query_selection.tsv",
        [{"query_id": k, "reason": v} for k, v in chosen.items()],
        ["query_id", "reason"],
    )
    markers = table(source / "archaeal_marker_hits.tsv")
    write_table(
        root / "archaeal_marker_hits.tsv",
        [r for r in markers if r["query_id"] in chosen],
        list(markers[0]),
    )
    (root / "selection.json").write_text(
        json.dumps(
            {
                "genes": len(genes),
                "bins": len({r["bin_id"] for r in genes}),
                "flagged_contigs": len(flagged),
                "screen_genes_sha256": sha256(source / "genes.faa"),
                "selection_criteria": config["selection_criteria"],
            },
            indent=2,
        )
        + "\n"
    )
    print((root / "selection.json").read_text())


if __name__ == "__main__":
    main()
