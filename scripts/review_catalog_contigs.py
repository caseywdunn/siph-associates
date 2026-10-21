#!/usr/bin/env python3
"""Record domain/contig decisions from completed protein and genome comparisons."""

import json
from collections import Counter, defaultdict
from pathlib import Path

from assembly_completion import table, unique, write_table

ROOT = Path("data/results/post_expansion/2026-10-20")


def main() -> None:
    bins = unique(table("figures/catalog_qc/domain_review.tsv"), "mag_id")
    contigs = unique(table(ROOT / "qc/contig_evidence.tsv"), "contig_id")
    nr = defaultdict(list)
    for row in table(
        "data/results/catalog_qc/2026-10-20/nr_followup/review/genes_nr.tsv"
    ):
        nr[row["contig_id"]].append(row)
    alignments = defaultdict(list)
    for ref in ("P_physalis", "N_septata", "H_rubrum", "Perkinsela"):
        with (ROOT / f"qc/curation/candidate_vs_{ref}.paf").open() as handle:
            for line in handle:
                x = line.rstrip().split("\t")
                alignments[x[0]].append(
                    (ref, (int(x[3]) - int(x[2])) / int(x[1]), int(x[9]) / int(x[10]))
                )
    decisions = []
    for key, row in contigs.items():
        b = bins[row["mag_id"]]
        if int(b["genes_Bacteria"]) <= int(b["genes_Eukaryota"]):
            continue
        counts = Counter(r["assignment"] for r in nr[key])
        euk = [r for r in nr[key] if r["assignment"] == "Eukaryota"]
        flagged = int(row["genes_Eukaryota"]) >= 3 and int(
            row["genes_Eukaryota"]
        ) > int(row["genes_Bacteria"])
        if not (euk or counts["Archaea"] or flagged or alignments[key]):
            continue
        # These comparisons concern assembled fragments, not read presence.
        matches = [x for x in alignments[key] if x[1] >= 0.5 and x[2] >= 0.75]
        if (
            matches
            and int(row["genes_Bacteria"]) <= 1
            and int(row["tiles_Bacteria"]) <= 1
            and counts["Bacteria"] == 0
        ):
            decision = "remove"
            reason = "Broad nucleotide match to a host reference or excluded eukaryotic bin, without coherent bacterial protein support."
        elif (
            euk
            and int(row["genes_Bacteria"]) == 0
            and int(row["tiles_Bacteria"]) == 0
            and counts["Bacteria"] == 0
            and (len(euk) >= 2 or max(int(r["best_alignment_aa"]) for r in euk) >= 150)
        ):
            decision = "remove"
            reason = "Eukaryotic protein and translated-window evidence; nr confirms a substantial protein match, with no bacterial assignment on this fragment."
        else:
            decision = "retain"
            reason = "Bacterial context or nr reassignment supports retention; isolated short/conserved cross-domain matches do not establish foreign contig origin."
        decisions.append(
            {
                "mag_id": row["mag_id"],
                "contig_id": key,
                "length": row["length"],
                "source_bin_sha256": b["sha256"],
                "decision": decision,
                "reason": reason,
                "swissprot_bacterial_genes": row["genes_Bacteria"],
                "swissprot_eukaryotic_genes": row["genes_Eukaryota"],
                "nr_bacterial_genes": counts["Bacteria"],
                "nr_eukaryotic_genes": counts["Eukaryota"],
                "nr_archaeal_genes": counts["Archaea"],
                "nr_eukaryotic_hits": "; ".join(r["best_title"] for r in euk),
                "genomic_matches": "; ".join(
                    f"{ref}:fraction={fraction:.4f},identity={identity:.4f}"
                    for ref, fraction, identity in alignments[key]
                ),
            }
        )
    write_table(
        "figures/catalog_qc/contig_decisions.tsv", decisions, list(decisions[0])
    )
    bin_decisions = []
    for key, row in bins.items():
        bacterial = int(row["genes_Bacteria"]) > int(row["genes_Eukaryota"])
        removed = [
            r for r in decisions if r["mag_id"] == key and r["decision"] == "remove"
        ]
        bin_decisions.append(
            {
                "mag_id": key,
                "source_sha256": row["sha256"],
                "decision": "retain_subject_to_quality" if bacterial else "exclude",
                "reason": "Broad bacterial evidence and GTDB markers; check quality after contig curation."
                if bacterial
                else "Eukaryote-enriched evidence without a supported prokaryotic genome assignment; includes the four known false-archaeal controls.",
                "removed_contigs": len(removed),
                "removed_bases": sum(int(r["length"]) for r in removed),
            }
        )
    write_table(
        "figures/catalog_qc/bin_decisions.tsv", bin_decisions, list(bin_decisions[0])
    )
    print(
        json.dumps(
            {
                "reviewed_contigs": len(decisions),
                "contig_decisions": dict(Counter(r["decision"] for r in decisions)),
                "bin_decisions": dict(Counter(r["decision"] for r in bin_decisions)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
