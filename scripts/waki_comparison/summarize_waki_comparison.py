#!/usr/bin/env python3
"""Make a specimen/marker comparison from retained BLAST alignments.

Every specimen and marker has a row. A top local hit is not a taxonomic call.
Different haplotypes of each reference taxon are compared on the same subject
coordinates before selecting the best competing reference.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from compare_waki_assemblies import read_tsv, sha256, subject_map, write_tsv


def comparison_taxon(row: dict[str, str]) -> str:
    """Retain Waki's grouped taxa and its explicit earlier-record identifications."""
    name = row["paper_taxon"].replace("_", " ").strip()
    if name.startswith("Hemiuridae sp. A"):
        return "Dinurus barbatus"
    if name.startswith("Sclerodistomidae sp. A"):
        return "Bathycotyle branchialis"
    return name


def shared_comparison(hits: list[dict]) -> tuple[list[dict], list[int]]:
    """Score every haplotype on identical forward subject-base coordinates.

    Subject bases represented by an insertion relative to the reference count
    against identity. Reference insertions have no subject coordinate and are
    excluded; full BLAST gap-inclusive identities remain separately reported.
    """
    maps = [subject_map(hit) for hit in hits]
    positions = sorted(set.intersection(*(set(mapped) for mapped in maps)))
    positions = [p for p in positions if all(mapped[p][1] in "ACGT" for mapped in maps)]
    scores = []
    for hit, mapped in zip(hits, maps, strict=True):
        matches = sum(mapped[p][0] == mapped[p][1] for p in positions)
        scores.append(
            {
                "hit": hit,
                "matches": matches,
                "positions": len(positions),
                "identity": 100 * matches / len(positions) if positions else None,
            }
        )
    return scores, positions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--analysis", type=Path, default=Path("data/results/waki_comparison")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("data/results/waki_comparison/report")
    )
    args = parser.parse_args()
    out = args.output
    if out.exists() and any(out.iterdir()):
        raise SystemExit(
            "Choose a fresh report output directory to preserve previous results"
        )
    out.mkdir(parents=True, exist_ok=True)
    receipt_path = args.analysis / "provenance.json"
    receipt = json.loads(receipt_path.read_text())
    if receipt["status"] != "complete":
        raise ValueError("Underlying BLAST analysis is incomplete")
    config = receipt["resolved_parameters"]
    inputs = [
        receipt_path,
        args.analysis / "reference_assembly_hits.tsv",
        args.analysis / "assemblies.fasta",
        Path(config["reference_metadata"]),
        Path(config["context_metadata"]),
        Path(config["manifest"]),
    ]
    for path in inputs[1:3]:
        if receipt["outputs_sha256"][str(path)] != sha256(path):
            raise ValueError(f"Analysis output checksum changed: {path}")
    metadata = {}
    for path in inputs[3:5]:
        for row in read_tsv(path):
            metadata[row["accession_version"]] = row
    with inputs[5].open() as handle:
        host_metadata = {
            row["library_id"].replace(":", "__"): row["species_current"]
            for row in csv.DictReader(handle)
        }
    all_hits = read_tsv(inputs[1])
    by_specimen_marker = defaultdict(list)
    for hit in all_hits:
        hit["comparison_taxon"] = comparison_taxon(metadata[hit["qseqid"]])
        by_specimen_marker[(hit["sample_id"], hit["marker"])].append(hit)
    rows = []
    detailed = []
    segments = []
    aligned = []
    candidate_contigs = set()
    assemblies = {record.id: record for record in SeqIO.parse(inputs[2], "fasta")}
    for sample in config["samples"]:
        for marker in ("28S", "ITS2", "COI"):
            hits = by_specimen_marker[(sample, marker)]
            waki = [h for h in hits if h["reference_set"] == "Waki2026"]
            row = {
                "sample_id": sample,
                "host": host_metadata[sample],
                "marker": marker,
                "result": "no_significant_alignment_to_Waki_records",
            }
            if not waki:
                rows.append(row)
                continue
            focal = max(waki, key=lambda h: float(h["bitscore"]))
            ref = metadata[focal["qseqid"]]
            mapped = subject_map(focal)
            eligible = [
                h
                for h in hits
                if h["sseqid"] == focal["sseqid"]
                and len(set(subject_map(h)) & set(mapped)) / len(mapped)
                >= config["competitor_min_locus_overlap_fraction"]
            ]
            scores, positions = shared_comparison(eligible)
            if not positions:
                raise ValueError(
                    f"No shared positions in comparable local alignments: {sample}, {marker}"
                )
            focal_score = next(
                s for s in scores if s["hit"]["hit_id"] == focal["hit_id"]
            )
            same_taxon = [
                s
                for s in scores
                if s["hit"]["comparison_taxon"] == focal["comparison_taxon"]
            ]
            other_taxa = [
                s
                for s in scores
                if s["hit"]["comparison_taxon"] != focal["comparison_taxon"]
            ]
            best_same = max(
                same_taxon, key=lambda s: (s["identity"], float(s["hit"]["bitscore"]))
            )
            competitor = (
                max(
                    other_taxa,
                    key=lambda s: (s["identity"], float(s["hit"]["bitscore"])),
                )
                if other_taxa
                else None
            )
            long = int(focal["length"]) >= config["long_homology_min_alignment_bp"]
            row.update(
                {
                    "host": focal["host"],
                    "result": "long_marker_homology"
                    if long
                    else "short_local_alignment_only",
                    "waki_accession": focal["qseqid"],
                    "waki_taxon": focal["comparison_taxon"],
                    "waki_genbank_organism": focal["reference_taxon"],
                    "waki_host": ref["host"],
                    "waki_developmental_stage": ref["developmental_stage"],
                    "assembly_contig": focal["sseqid"],
                    "assembly_contig_length_bp": focal["slen"],
                    "reference_length_bp": focal["qlen"],
                    "alignment_length_bp": focal["length"],
                    "identical_bases": sum(
                        q == s and q != "-"
                        for q, s in zip(focal["qseq"], focal["sseq"], strict=True)
                    ),
                    "alignment_identity_percent": 100
                    * sum(
                        q == s and q != "-"
                        for q, s in zip(focal["qseq"], focal["sseq"], strict=True)
                    )
                    / int(focal["length"]),
                    "reference_coverage_percent": focal["query_coverage_percent"],
                    "reference_start": focal["qstart"],
                    "reference_end": focal["qend"],
                    "subject_start": focal["sstart"],
                    "subject_end": focal["send"],
                    "bitscore": focal["bitscore"],
                    "evalue": focal["evalue"],
                    "physical_link_to_accepted_18S": "none_found"
                    if not focal["same_specimen_ssu_sequences"]
                    else focal["same_specimen_ssu_sequences"],
                    "shared_subject_start": min(positions),
                    "shared_subject_end": max(positions),
                    "shared_subject_bases": len(positions),
                    "focal_identity_shared_percent": round(focal_score["identity"], 4),
                    "best_same_taxon_accession_shared": best_same["hit"]["qseqid"],
                    "best_same_taxon_identity_shared_percent": round(
                        best_same["identity"], 4
                    ),
                    "eligible_comparison_taxa": len(
                        {h["comparison_taxon"] for h in eligible}
                    ),
                    "interpretation": "reference_affinity_only; not_a_species_assignment",
                    "marker_caveat": "5.8S/ITS2/28S amplicon; internal boundaries unannotated"
                    if marker == "ITS2"
                    else "",
                }
            )
            if competitor:
                row.update(
                    {
                        "competing_taxon": competitor["hit"]["comparison_taxon"],
                        "competing_accession": competitor["hit"]["qseqid"],
                        "competing_reference_set": competitor["hit"]["reference_set"],
                        "competing_identity_shared_percent": round(
                            competitor["identity"], 4
                        ),
                    }
                )
            for score in scores:
                hit = score["hit"]
                detailed.append(
                    {
                        "sample_id": sample,
                        "marker": marker,
                        "focal_hit_id": focal["hit_id"],
                        "hit_id": hit["hit_id"],
                        "accession": hit["qseqid"],
                        "reference_set": hit["reference_set"],
                        "comparison_taxon": hit["comparison_taxon"],
                        "genbank_organism": hit["reference_taxon"],
                        "subject_id": hit["sseqid"],
                        "shared_subject_start": min(positions),
                        "shared_subject_end": max(positions),
                        "shared_subject_bases": len(positions),
                        "identical_bases": score["matches"],
                        "identity_shared_percent": round(score["identity"], 4),
                        "blast_identity_percent": hit["pident"],
                        "blast_alignment_length": hit["length"],
                    }
                )
            if long:
                candidate_contigs.add(focal["sseqid"])
                lo, hi = sorted((int(focal["sstart"]), int(focal["send"])))
                sequence = assemblies[focal["sseqid"]].seq[lo - 1 : hi]
                if int(focal["send"]) < int(focal["sstart"]):
                    sequence = sequence.reverse_complement()
                segment_id = sample + "_" + marker + "_" + focal["hit_id"]
                segments.append(
                    SeqRecord(
                        sequence,
                        id=segment_id,
                        description=f"{focal['sseqid']}:{focal['sstart']}-{focal['send']} oriented_to_{focal['qseqid']}",
                    )
                )
                aligned.extend(
                    [
                        SeqRecord(
                            Seq(focal["qseq"]),
                            id=segment_id + "_reference",
                            description=focal["qseqid"],
                        ),
                        SeqRecord(
                            Seq(focal["sseq"]),
                            id=segment_id + "_assembly",
                            description=focal["sseqid"],
                        ),
                    ]
                )
            rows.append(row)
    fields = [
        "sample_id",
        "host",
        "marker",
        "result",
        "waki_accession",
        "waki_taxon",
        "waki_genbank_organism",
        "waki_host",
        "waki_developmental_stage",
        "assembly_contig",
        "assembly_contig_length_bp",
        "reference_length_bp",
        "alignment_length_bp",
        "identical_bases",
        "alignment_identity_percent",
        "reference_coverage_percent",
        "reference_start",
        "reference_end",
        "subject_start",
        "subject_end",
        "bitscore",
        "evalue",
        "physical_link_to_accepted_18S",
        "shared_subject_start",
        "shared_subject_end",
        "shared_subject_bases",
        "focal_identity_shared_percent",
        "best_same_taxon_accession_shared",
        "best_same_taxon_identity_shared_percent",
        "competing_taxon",
        "competing_accession",
        "competing_reference_set",
        "competing_identity_shared_percent",
        "eligible_comparison_taxa",
        "interpretation",
        "marker_caveat",
    ]
    write_tsv(out / "specimen_marker_summary.tsv", rows, fields)
    detail_fields = [
        "sample_id",
        "marker",
        "focal_hit_id",
        "hit_id",
        "accession",
        "reference_set",
        "comparison_taxon",
        "genbank_organism",
        "subject_id",
        "shared_subject_start",
        "shared_subject_end",
        "shared_subject_bases",
        "identical_bases",
        "identity_shared_percent",
        "blast_identity_percent",
        "blast_alignment_length",
    ]
    write_tsv(
        out / "all_haplotype_shared_interval_comparisons.tsv", detailed, detail_fields
    )
    SeqIO.write(segments, out / "top_marker_fragments.fasta", "fasta")
    SeqIO.write(aligned, out / "top_marker_pairwise_alignments.fasta", "fasta")
    SeqIO.write(
        (assemblies[key] for key in sorted(candidate_contigs)),
        out / "top_marker_contigs.fasta",
        "fasta",
    )
    source = Path(__file__)
    helper = source.with_name("compare_waki_assemblies.py")
    shutil.copyfile(source, out / source.name)
    record = {
        "status": "complete",
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "command": [sys.executable] + sys.argv,
        "source_snapshot": str(out / source.name),
        "script_sha256": sha256(source),
        "helper_sha256": {str(helper): sha256(helper)},
        "inputs_sha256": {str(path): sha256(path) for path in inputs},
        "parameters": {
            "min_competitor_overlap_fraction": config[
                "competitor_min_locus_overlap_fraction"
            ],
            "same_interval_comparison": "all eligible reference haplotypes scored before choosing best per taxon",
            "identity_denominator": "same unambiguous subject coordinates, excluding reference insertions",
            "taxon_aliases_from_Waki": {
                "Hemiuridae sp. A": "Dinurus barbatus",
                "Sclerodistomidae sp. A": "Bathycotyle branchialis",
            },
        },
        "checks": {
            "eight_specimens_three_markers": len(rows) == 24,
            "number_of_fragment_alignments": len(segments),
        },
    }
    record["outputs_sha256"] = {
        str(path): sha256(path) for path in out.iterdir() if path.is_file()
    }
    (out / "provenance.json").write_text(json.dumps(record, indent=2) + "\n")
    print(
        f"Wrote {len(rows)} specimen/marker rows and {len(detailed)} shared-interval haplotype comparisons"
    )


if __name__ == "__main__":
    main()
