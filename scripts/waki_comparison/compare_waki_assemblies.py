#!/usr/bin/env python3
"""Compare Waki's deposited markers with existing siphonophore assemblies.

BLAST retains all significant local alignments, including mismatches to the
SSU marker. Long matches are homology candidates, never species calls. The
comparison changes no accepted presence grades.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

BLAST_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
    "qseq",
    "sseq",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_fasta(path: Path) -> dict[str, SeqRecord]:
    records = {}
    for record in SeqIO.parse(path, "fasta"):
        if not record.seq or record.id in records:
            raise ValueError(
                f"Empty sequence or duplicate identifier in {path}: {record.id}"
            )
        if set(str(record.seq).upper()) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Unexpected nucleotide alphabet in {path}: {record.id}")
        records[record.id] = record
    return records


def subject_map(hit: dict) -> dict[int, tuple[str, str]]:
    """Map 1-based subject coordinates to paired bases in forward orientation."""
    position = int(hit["sstart"])
    direction = 1 if int(hit["send"]) >= position else -1
    mapped = {}
    for query, subject in zip(hit["qseq"], hit["sseq"], strict=True):
        if subject != "-":
            if direction == -1:
                query, subject = (
                    str(Seq(query).complement()),
                    str(Seq(subject).complement()),
                )
            mapped[position] = (query.upper(), subject.upper())
            position += direction
    if position - direction != int(hit["send"]):
        raise ValueError("BLAST alignment sequence and subject coordinates disagree")
    return mapped


def read_blast(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, fieldnames=BLAST_FIELDS, delimiter="\t"))
    for row in rows:
        if len(row["qseq"]) != int(row["length"]) or len(row["sseq"]) != int(
            row["length"]
        ):
            raise ValueError(f"BLAST aligned length mismatch in {path}")
        subject_map(row)
    return rows


def marker_from_record(record: SeqRecord) -> str:
    description = record.description.lower()
    if "its2" in description or "internal transcribed spacer" in description:
        return "ITS2"
    if "28s" in description:
        return "28S"
    if "cytochrome" in description or "cox1" in description or "coi" in description:
        return "COI"
    raise ValueError(
        f"Cannot recognize reference marker: {record.id}: {record.description}"
    )


def summarize(
    out: Path,
    config: dict,
    reference_info: dict,
    subjects: dict,
    ssu_info: dict,
    manifest: dict,
) -> dict:
    assembly_hits = read_blast(out / "references_vs_assemblies.blast.tsv")
    ssu_hits = read_blast(out / "accepted_ssu_vs_assemblies.blast.tsv")
    direct_hits = read_blast(out / "references_vs_accepted_ssu.blast.tsv")
    links = []
    by_contig = defaultdict(list)
    for hit in ssu_hits:
        query = ssu_info[hit["qseqid"]]
        target = subjects[hit["sseqid"]]
        coverage = (
            100 * (abs(int(hit["qend"]) - int(hit["qstart"])) + 1) / int(hit["qlen"])
        )
        linked = (
            query["sample_id"] == target["sample_id"]
            and float(hit["pident"]) >= config["physical_ssu_link_min_identity_percent"]
            and coverage >= config["physical_ssu_link_min_query_coverage_percent"]
        )
        row = {
            **hit,
            "source_sequence_id": query["sequence_id"],
            "source_sample_id": query["sample_id"],
            "target_sample_id": target["sample_id"],
            "query_coverage_percent": round(coverage, 3),
            "same_specimen_physical_link": str(linked).lower(),
        }
        links.append(row)
        if linked:
            by_contig[hit["sseqid"]].append(query["sequence_id"])
    write_tsv(
        out / "ssu_assembly_links.tsv",
        links,
        BLAST_FIELDS
        + [
            "source_sequence_id",
            "source_sample_id",
            "target_sample_id",
            "query_coverage_percent",
            "same_specimen_physical_link",
        ],
    )

    annotated = []
    for index, hit in enumerate(assembly_hits, 1):
        reference = reference_info[hit["qseqid"]]
        sample = subjects[hit["sseqid"]]["sample_id"]
        query_coverage = (
            100 * (abs(int(hit["qend"]) - int(hit["qstart"])) + 1) / int(hit["qlen"])
        )
        long = int(hit["length"]) >= config["long_homology_min_alignment_bp"]
        row = {
            "hit_id": f"H{index:06d}",
            **hit,
            "sample_id": sample,
            "host": manifest[sample]["species_current"],
            "reference_taxon": reference["taxon"],
            "marker": reference["marker"],
            "reference_set": reference["reference_set"],
            "query_coverage_percent": round(query_coverage, 3),
            "long_homology_candidate": str(long).lower(),
            "same_specimen_ssu_sequences": ";".join(by_contig[hit["sseqid"]]),
            "marker_region": "unpartitioned_5.8S_ITS2_28S_amplicon"
            if reference["marker"] == "ITS2"
            else reference["marker"],
        }
        annotated.append(row)
    extra = [
        "sample_id",
        "host",
        "reference_taxon",
        "marker",
        "reference_set",
        "query_coverage_percent",
        "long_homology_candidate",
        "same_specimen_ssu_sequences",
        "marker_region",
    ]
    write_tsv(
        out / "reference_assembly_hits.tsv",
        annotated,
        ["hit_id"] + BLAST_FIELDS + extra,
    )
    long_hits = [row for row in annotated if row["long_homology_candidate"] == "true"]
    write_tsv(
        out / "long_homology_candidates.tsv",
        long_hits,
        ["hit_id"] + BLAST_FIELDS + extra,
    )

    # Every eligible competitor must cover at least 80% of the focal subject
    # interval. Compare all of them over exactly their shared subject bases.
    comparison = []
    loci = defaultdict(list)
    for row in long_hits:
        loci[(row["sseqid"], row["marker"])].append(row)
    for (contig, marker), hits in sorted(loci.items()):
        remaining = sorted(hits, key=lambda h: -float(h["bitscore"]))
        locus_n = 0
        while remaining:
            locus_n += 1
            focal = remaining[0]
            focal_map = subject_map(focal)
            comparable = []
            remainder = []
            for hit in remaining:
                overlap = len(set(subject_map(hit)) & set(focal_map)) / len(focal_map)
                (
                    comparable
                    if overlap >= config["competitor_min_locus_overlap_fraction"]
                    else remainder
                ).append(hit)
            remaining = remainder
            best_by_taxon = {}
            for hit in comparable:
                taxon = hit["reference_taxon"]
                if taxon not in best_by_taxon:
                    best_by_taxon[taxon] = hit
            maps = [subject_map(h) for h in best_by_taxon.values()]
            common = set.intersection(*(set(m) for m in maps))
            if not common:
                continue
            for hit in best_by_taxon.values():
                mapped = subject_map(hit)
                eligible = [p for p in common if mapped[p][1] in "ACGT"]
                matches = sum(mapped[p][0] == mapped[p][1] for p in eligible)
                ambiguous = sum(mapped[p][0] not in "ACGT-" for p in eligible)
                comparison.append(
                    {
                        "locus_id": f"{contig}:{marker}:{locus_n}",
                        "sample_id": hit["sample_id"],
                        "subject_id": contig,
                        "marker": marker,
                        "reference_taxon": hit["reference_taxon"],
                        "reference_accession": hit["qseqid"],
                        "reference_set": hit["reference_set"],
                        "hit_id": hit["hit_id"],
                        "shared_subject_start": min(common),
                        "shared_subject_end": max(common),
                        "shared_subject_bases": len(common),
                        "compared_subject_ACGT_bases": len(eligible),
                        "reference_ambiguous_bases": ambiguous,
                        "identical_bases": matches,
                        "identity_on_shared_subject_bases_percent": round(
                            100 * matches / len(eligible), 4
                        )
                        if eligible
                        else "",
                        "competing_reference_taxa": len(best_by_taxon),
                        "limitation": "subject-coordinate comparison excludes reference insertions; ITS2 boundaries unannotated"
                        if marker == "ITS2"
                        else "subject-coordinate comparison excludes reference insertions",
                    }
                )
    comparison_fields = [
        "locus_id",
        "sample_id",
        "subject_id",
        "marker",
        "reference_taxon",
        "reference_accession",
        "reference_set",
        "hit_id",
        "shared_subject_start",
        "shared_subject_end",
        "shared_subject_bases",
        "compared_subject_ACGT_bases",
        "reference_ambiguous_bases",
        "identical_bases",
        "identity_on_shared_subject_bases_percent",
        "competing_reference_taxa",
        "limitation",
    ]
    write_tsv(out / "same_interval_comparisons.tsv", comparison, comparison_fields)
    summary = {
        "assemblies": len(config["samples"]),
        "reference_sequences": len(reference_info),
        "reference_counts": dict(
            Counter(
                (r["reference_set"] + ":" + r["marker"])
                for r in reference_info.values()
            )
        ),
        "accepted_helminth_ssu_sequences": len(ssu_info),
        "assembly_blast_hsps": len(assembly_hits),
        "reference_vs_accepted_ssu_hsps": len(direct_hits),
        "ssu_vs_assembly_hsps": len(ssu_hits),
        "long_homology_candidates": len(long_hits),
        "long_homology_samples": sorted({h["sample_id"] for h in long_hits}),
        "long_homology_marker_counts": dict(Counter(h["marker"] for h in long_hits)),
        "physically_linked_ssu_candidates": sum(
            bool(h["same_specimen_ssu_sequences"]) for h in long_hits
        ),
        "same_interval_loci": len({r["locus_id"] for r in comparison}),
        "interpretation": "Local sequence homology against a selected reference panel. Not species assignments or new presence grades.",
        "ITS2_limit": "GenBank labels entire 5.8S/ITS2/28S amplicons without internal boundaries; identity cannot be assigned specifically to variable ITS2.",
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/waki_comparison.json")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("data/results/waki_comparison")
    )
    parser.add_argument("--include-context", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    if (out / "provenance.json").exists():
        raise SystemExit(
            "Choose a fresh --output directory to preserve prior execution records"
        )
    start = datetime.now(timezone.utc).isoformat()
    inputs = [args.config] + [
        Path(config[k])
        for k in (
            "references",
            "reference_genbank",
            "reference_metadata",
            "manifest",
            "grades",
            "assembled_ssu_table",
            "assembled_ssu_fasta",
        )
    ]
    reference_sets = [
        ("Waki2026", "references", "reference_genbank", "reference_metadata")
    ]
    if args.include_context:
        reference_sets.append(
            ("context", "context_references", "context_genbank", "context_metadata")
        )
        inputs.extend(
            Path(config[k])
            for k in ("context_references", "context_genbank", "context_metadata")
        )
    with Path(config["manifest"]).open() as handle:
        metadata = {
            r["library_id"].replace(":", "__"): r for r in csv.DictReader(handle)
        }
    refs = {}
    ref_info = {}
    features = []
    for label, fasta_key, gb_key, metadata_key in reference_sets:
        fasta = load_fasta(Path(config[fasta_key]))
        marker_lookup = {
            r["accession_version"]: r["marker"]
            for r in read_tsv(Path(config[metadata_key]))
        }
        for record in SeqIO.parse(config[gb_key], "genbank"):
            if (
                record.id not in fasta
                or str(record.seq).upper() != str(fasta[record.id].seq).upper()
            ):
                raise ValueError(f"GenBank/FASTA disagreement: {record.id}")
            if record.id in refs:
                if str(record.seq).upper() != str(refs[record.id].seq).upper():
                    raise ValueError(f"Reference sets disagree: {record.id}")
                continue
            refs[record.id] = record
            ref_info[record.id] = {
                "taxon": record.annotations["organism"],
                "marker": marker_lookup[record.id],
                "reference_set": label,
            }
            for feature in record.features:
                features.append(
                    {
                        "accession_version": record.id,
                        "feature_type": feature.type,
                        "start_0based": int(feature.location.start),
                        "end_exclusive": int(feature.location.end),
                        "location": str(feature.location),
                        "qualifiers": json.dumps(feature.qualifiers, sort_keys=True),
                    }
                )
    SeqIO.write(refs.values(), out / "search_references.fasta", "fasta")
    write_tsv(
        out / "reference_features.tsv",
        features,
        [
            "accession_version",
            "feature_type",
            "start_0based",
            "end_exclusive",
            "location",
            "qualifiers",
        ],
    )
    grades = read_tsv(Path(config["grades"]))
    helminth_libraries = {
        r["sample_id"]
        for r in grades
        if r["count_as_detection"] == "true"
        and r["role"] == "parasite"
        and r["reporting_unit"] == "Platyhelminthes"
    }
    assemblies = read_tsv(Path(config["assembled_ssu_table"]))
    ssu_fasta = load_fasta(Path(config["assembled_ssu_fasta"]))
    ssu_info = {}
    ssu_records = []
    for row in assemblies:
        if (
            row["sample_id"] in helminth_libraries
            and "Platyhelminthes" in row["lineage"]
        ):
            alias = f"SSU{len(ssu_info) + 1:04d}"
            ssu_info[alias] = row
            ssu_records.append(
                SeqRecord(ssu_fasta[row["sequence_id"]].seq, id=alias, description="")
            )
    SeqIO.write(ssu_records, out / "accepted_helminth_ssu.fasta", "fasta")
    write_tsv(
        out / "ssu_identifiers.tsv",
        [{"alias": a, **r} for a, r in ssu_info.items()],
        ["alias"] + list(next(iter(ssu_info.values()))),
    )
    subjects = {}
    inventory = []
    with (out / "assemblies.fasta").open("w") as combined:
        for sample in config["samples"]:
            path = Path(config["assembly_pattern"].format(sample_id=sample))
            inputs.append(path)
            records = load_fasta(path)
            for record in records.values():
                if record.id in subjects or not record.id.startswith(sample + "__"):
                    raise ValueError(
                        f"Assembly identifier is not globally unique: {record.id}"
                    )
                subjects[record.id] = {"sample_id": sample, "length": len(record)}
            SeqIO.write(records.values(), combined, "fasta")
            inventory.append(
                {
                    "sample_id": sample,
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "contigs": len(records),
                    "nucleotides": sum(len(r) for r in records.values()),
                }
            )
    write_tsv(
        out / "assembly_inventory.tsv",
        inventory,
        ["sample_id", "path", "size_bytes", "contigs", "nucleotides"],
    )
    script = Path(__file__)
    shutil.copyfile(script, out / "compare_waki_assemblies.py")
    shutil.copyfile(args.config, out / "config.json")
    receipt = {
        "status": "running",
        "started_utc": start,
        "command": [sys.executable] + sys.argv,
        "run_id": start,
        "inputs_sha256": {str(p): sha256(p) for p in inputs},
        "source_snapshot": str(out / script.name),
        "script_sha256": sha256(script),
        "python": sys.version,
        "biopython": Bio.__version__,
        "platform": platform.platform(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "resolved_parameters": config,
        "tools": {},
        "commands": [],
    }
    for tool in ("blastn", "makeblastdb"):
        receipt["tools"][tool] = subprocess.check_output(
            [config[tool], "-version"], text=True
        ).strip()
    (out / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")

    def run(command: list[str], label: str) -> None:
        print(label, flush=True)
        receipt["commands"].append(command)
        (out / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")
        with (out / f"{label}.log").open("w") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)

    try:
        run(
            [
                config["makeblastdb"],
                "-in",
                str(out / "assemblies.fasta"),
                "-dbtype",
                "nucl",
                "-out",
                str(out / "assembly_db"),
            ],
            "makeblastdb",
        )
        params = config["blast_parameters"]
        common = [
            "-task",
            params["task"],
            "-word_size",
            str(params["word_size"]),
            "-dust",
            params["dust"],
            "-evalue",
            params["evalue"],
            "-max_target_seqs",
            str(params["max_target_seqs"]),
            "-num_threads",
            str(params["threads"]),
            "-outfmt",
            "6 " + " ".join(BLAST_FIELDS),
        ]
        for query, target, output, label in [
            (
                out / "search_references.fasta",
                ["-db", str(out / "assembly_db")],
                "references_vs_assemblies.blast.tsv",
                "references_vs_assemblies",
            ),
            (
                out / "accepted_helminth_ssu.fasta",
                ["-db", str(out / "assembly_db")],
                "accepted_ssu_vs_assemblies.blast.tsv",
                "accepted_ssu_vs_assemblies",
            ),
            (
                out / "search_references.fasta",
                ["-subject", str(out / "accepted_helminth_ssu.fasta")],
                "references_vs_accepted_ssu.blast.tsv",
                "references_vs_accepted_ssu",
            ),
        ]:
            run(
                [config["blastn"], "-query", str(query)]
                + target
                + common
                + ["-out", str(out / output)],
                label,
            )
        summary = summarize(out, config, ref_info, subjects, ssu_info, metadata)
        receipt["summary"] = summary
        receipt["status"] = "complete"
        print(json.dumps(summary, indent=2), flush=True)
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error"] = repr(exc)
        raise
    finally:
        receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
        receipt["outputs_sha256"] = {
            str(p): sha256(p)
            for p in out.iterdir()
            if p.is_file() and p.name != "provenance.json"
        }
        (out / "provenance.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    main()
