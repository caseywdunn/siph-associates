#!/usr/bin/env python3
"""Apply recorded bin/contig decisions and select rechecked viral regions."""

import argparse
import json
import shutil
from pathlib import Path

from assembly_completion import sha256, table, unique, write_table
from Bio import SeqIO

ROOT = Path("data/results/post_expansion/2026-10-20")
OUT = ROOT / "curated"


def prepare() -> None:
    decisions = unique(table("figures/catalog_qc/bin_decisions.tsv"), "mag_id")
    contigs = unique(table("figures/catalog_qc/contig_decisions.tsv"), "contig_id")
    candidates = unique(table(ROOT / "qc/candidates.tsv"), "mag_id")
    if set(candidates) != set(decisions):
        raise ValueError("Candidate and decision identities disagree")
    dest = OUT / "phase3_catalog/mags"
    (dest / "candidates").mkdir(parents=True, exist_ok=True)
    (dest / "changed").mkdir(exist_ok=True)
    records = []
    for key, d in decisions.items():
        source = ROOT / "qc/candidates" / f"{key}.fa"
        if sha256(source) != d["source_sha256"]:
            raise ValueError(f"Reviewed bin changed: {key}")
        if d["decision"] == "exclude":
            continue
        if d["decision"] != "retain_subject_to_quality":
            raise ValueError(f"Unreviewed bin: {key}")
        excluded = {
            r["contig_id"]
            for r in contigs.values()
            if r["mag_id"] == key and r["decision"] == "remove"
        }
        genome = list(SeqIO.parse(source, "fasta"))
        if not excluded <= {r.id for r in genome}:
            raise ValueError(f"Missing reviewed contig in {key}")
        target = dest / "candidates" / f"{key}.fa"
        kept = [r for r in genome if r.id not in excluded]
        if not kept:
            raise ValueError(f"Curation removed all of {key}")
        if excluded:
            SeqIO.write(kept, target, "fasta")
        else:
            shutil.copyfile(source, target)
        if excluded:
            changed = dest / "changed" / target.name
            changed.symlink_to(target.resolve())
        records.append(
            dict(
                candidates[key],
                curated_sha256=sha256(target),
                removed_contigs=len(excluded),
                removed_bases=sum(len(r) for r in genome if r.id in excluded),
            )
        )
    write_table(dest / "before_quality.tsv", records, list(records[0]))
    # Old placements are reused only for the source-identical bins verified in staging.
    previous = unique(table("figures/catalog_qc/candidate_review.tsv"), "mag_id")
    new = table(ROOT / "qc/new_gtdb/gtdbtk.bac120.summary.tsv")
    placed = {}
    for kind in ("bac120", "ar53"):
        for r in table(f"data/results/phase3_catalog/mags/gtdbtk.{kind}.summary.tsv"):
            if (
                r["user_genome"] in previous
                and previous[r["user_genome"]]["previous_sequence_matches"] == "true"
            ):
                placed[r["user_genome"]] = r
    placed.update({r["user_genome"]: r for r in new})
    for kind, domain in (("bac120", "d__Bacteria"), ("ar53", "d__Archaea")):
        rows = [
            placed[r["mag_id"]]
            for r in records
            if placed[r["mag_id"]]["classification"].startswith(domain)
        ]
        write_table(dest / f"gtdbtk.{kind}.summary.tsv", rows, list(new[0]))
    print(
        f"Curated {len(records)} bacterial bins; {sum(bool(r['removed_contigs']) for r in records)} need updated CheckM2"
    )


def quality() -> None:
    dest = OUT / "phase3_catalog/mags"
    q = unique(table(dest / "changed_checkm2/quality_report.tsv"), "Name")
    policy = json.loads(Path("config/phase3_catalog.json").read_text())["mags"]
    rows = table(dest / "before_quality.tsv")
    changed = {r["mag_id"] for r in rows if int(r["removed_contigs"])}
    if set(q) != changed:
        raise ValueError("CheckM2 results do not cover exactly the edited genomes")
    for row in rows:
        if row["mag_id"] in q:
            updated = q[row["mag_id"]]
            for dest_field, src_field in [
                ("completeness", "Completeness"),
                ("contamination", "Contamination"),
                ("genome_size", "Genome_Size"),
                ("gc_content", "GC_Content"),
                ("contig_n50", "Contig_N50"),
                ("checkm2_model", "Completeness_Model_Used"),
            ]:
                row[dest_field] = updated[src_field]
        comp, contam = float(row["completeness"]), float(row["contamination"])
        row["quality_score"] = round(comp - 5 * contam, 2)
        row["near_complete"] = str(
            comp >= policy["near_complete_min_completeness"]
            and contam < policy["near_complete_max_contamination_exclusive"]
        ).lower()
        row["post_curation_quality"] = (
            "pass"
            if comp >= policy["min_completeness"]
            and contam < policy["max_contamination_exclusive"]
            else "exclude_quality"
        )
    write_table("figures/catalog_qc/curated_mag_quality.tsv", rows, list(rows[0]))
    keep = [r for r in rows if r["post_curation_quality"] == "pass"]
    write_table(dest / "candidates.tsv", keep, list(rows[0]))
    print(f"{len(keep)} of {len(rows)} curated genomes pass quality thresholds")


def viruses() -> None:
    root = ROOT / "qc/curation"
    checks = unique(table(root / "viral_checkv/quality_summary.tsv"), "contig_id")
    boundaries = unique(table(root / "viral_checkv/contamination.tsv"), "contig_id")
    rows = table(root / "viral_regions.tsv")
    rules = json.loads(Path("config/phase3_catalog.json").read_text())["viruses"]
    keep = []
    for r in rows:
        q = checks[r["contig_id"]]
        if boundaries[r["contig_id"]]["provirus"] != "No":
            raise ValueError(f"Unresolved cellular flanks: {r['contig_id']}")
        comp = float(q["completeness"]) if q["completeness"] not in ("", "NA") else 0
        r.update(
            checkv_quality=q["checkv_quality"],
            checkv_completeness=q["completeness"],
            checkv_completeness_method=q["completeness_method"],
            checkv_viral_genes=q["viral_genes"],
            checkv_host_genes=q["host_genes"],
            checkv_warnings=q["warnings"],
        )
        r["decision"] = "retain"
        if int(r["genomad_hallmarks"]) < rules["min_genomad_hallmarks"]:
            r["decision"] = "exclude_no_hallmark_in_retained_region"
        elif not (
            int(r["length"]) >= rules["min_length"]
            or comp >= rules["min_checkv_completeness"]
        ):
            r["decision"] = "exclude_below_length_and_completeness"
        if r["decision"] == "retain":
            keep.append(r)
    dest = OUT / "phase3_catalog/viruses"
    dest.mkdir(parents=True, exist_ok=True)
    write_table("figures/catalog_qc/viral_region_decisions.tsv", rows, list(rows[0]))
    write_table(dest / "included.tsv", keep, list(rows[0]))
    ids = {r["contig_id"] for r in keep}
    SeqIO.write(
        (r for r in SeqIO.parse(root / "viral_regions.fna", "fasta") if r.id in ids),
        dest / "included.fna",
        "fasta",
    )
    print(
        f"Retained {len(keep)} viral regions after boundary, hallmark and quality review"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "quality", "viruses"))
    args = parser.parse_args()
    {"prepare": prepare, "quality": quality, "viruses": viruses}[args.stage]()
