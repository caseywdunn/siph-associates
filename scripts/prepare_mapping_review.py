#!/usr/bin/env python3
"""Stage one assembly per library and an explicitly unapproved catalog-QC bundle."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

from assembly_completion import sha256, table, unique, write_table

PRODUCTS = {
    "assemblies": [".fasta"],
    "markers": [".gff"],
    "viruses": [".fna", ".tsv"],
    "checkv": [".tsv"],
    "bins": [".tar.gz", ".depth.tsv"],
    "checkm2": [".tsv", ".status"],
}


def check_hold(path: str = "config/post_expansion_remapping.json") -> None:
    hold = json.loads(Path(path).read_text())
    if (
        hold["status"] != "hold"
        or hold["automatic_remapping"]
        or not hold["requires_explicit_mapping_authorization"]
    ):
        raise ValueError(
            "Preparation must preserve the explicit full-cohort mapping hold"
        )


def mapping_input(row: dict, settings: dict) -> dict:
    """Inventory usable read derivatives; never regenerate or decompress raw reads."""
    sample, root = row["sample_id"], Path(row["input_root"])
    reference = row["host_route"] != "none"
    if reference:
        read_root = root / "phase3_cohort/inputs"
        provenance = root / "provenance/phase3_cohort/inputs" / f"{sample}.json"
        pairs = int(row["retained_pairs"])
    else:
        provenance = root / "provenance/trim" / f"{sample}.json"
        pairs = json.loads(provenance.read_text())["pairs_after_fastp"]
        candidates = [
            root / "retained/trimmed",
            root / "phase5_mapping/trimmed",
            Path("data/results/fish_12s/full_reads/2026-10-18/trimmed"),
            Path("/vast/palmer/scratch/dunn/cwd7/siph_associates/trimmed"),
        ]
        read_root = next(
            (
                p
                for p in candidates
                if all((p / f"{sample}_R{m}.fastq.gz").is_file() for m in (1, 2))
            ),
            None,
        )
    paths = [read_root / f"{sample}_R{m}.fastq.gz" for m in (1, 2)] if read_root else []
    available = len(paths) == 2 and all(p.is_file() for p in paths)
    if reference and not available:
        raise FileNotFoundError(f"Retained host-depleted reads missing: {sample}")
    if available:
        # Compare sizes with the completed producer's record before declaring ready.
        actual = json.loads(provenance.read_text())
        if not reference and "fish_12s/full_reads" in str(read_root):
            provenance = read_root.parent / "trim" / f"{sample}.json"
            actual = json.loads(provenance.read_text())
        elif not reference and "phase5_mapping" in str(read_root):
            provenance = root / "provenance/phase5_mapping/trim" / f"{sample}.json"
            actual = json.loads(provenance.read_text())
        if actual["sample_id"] != sample or len(actual["outputs"]) != 2:
            raise ValueError(f"Mapping input provenance mismatch: {sample}")
        if [p.stat().st_size for p in paths] != [
            o["size_bytes"] for o in actual["outputs"]
        ]:
            raise ValueError(f"Mapping input sizes differ from provenance: {sample}")
        if not reference and actual["pairs_after_fastp"] != pairs:
            raise ValueError(
                f"Mapping input pair count differs from baseline: {sample}"
            )
    return {
        "sample_id": sample,
        "species_current": row["species_current"],
        "host_route": row["host_route"],
        "mapping_read_source": "host_depleted" if reference else "full_capped_trimmed",
        "read_pairs": pairs,
        "availability": "ready"
        if available
        else "regenerate_trimmed_reads_when_mapping_is_authorized",
        "r1": str(paths[0]) if available else "",
        "r2": str(paths[1]) if available else "",
        "baseline_provenance": str(provenance),
        "baseline_provenance_sha256": sha256(provenance),
        "raw_r1_paths": row["r1_paths"],
        "raw_r2_paths": row["r2_paths"],
        "raw_pairs": row["raw_pairs"],
        "cap_pairs": 200000000,
        "cap_seed": 20260829,
    }


def stage(settings: dict) -> None:
    check_hold()
    cohort = table(settings["cohort"])
    if len(unique(cohort, "sample_id")) != settings["expected_cohort_samples"]:
        raise ValueError("Combined cohort must contain exactly 213 unique libraries")
    root = Path(settings["preparation_root"])
    summaries = {}
    for source in {r["result_root"] for r in cohort}:
        done = Path(source) / "stages/phase3_assembly.done"
        status = dict(line.split("\t", 1) for line in done.read_text().splitlines())
        if status["status"] != "PASS":
            raise ValueError(f"Incomplete assembly source: {source}")
        summaries[source] = unique(
            table(Path(source) / "phase3_cohort/assembly_summary.tsv"), "sample_id"
        )
    combined, inventory, mapping = [], [], []
    for row in cohort:
        sample, source = row["sample_id"], Path(row["result_root"]) / "phase3_cohort"
        summary = summaries[row["result_root"]][sample]
        if summary["assembly_status"] not in ("assembled", "no_contigs"):
            raise ValueError(f"Assembly was not attempted: {sample}")
        if any(
            str(summary[k]) != str(row[k])
            for k in ("host_route", "strategy", "retained_pairs")
        ):
            raise ValueError(
                f"Selected assembly does not match its read input: {sample}"
            )
        fasta = source / "assemblies" / f"{sample}.fasta"
        checksum = sha256(fasta)
        if checksum != summary["assembly_sha256"]:
            raise ValueError(f"Assembly checksum mismatch: {sample}")
        for kind, suffixes in PRODUCTS.items():
            for suffix in suffixes:
                origin = (source / kind / f"{sample}{suffix}").resolve(strict=True)
                destination = root / "phase3_cohort" / kind / origin.name
                destination.parent.mkdir(parents=True, exist_ok=True)
                if destination.is_symlink():
                    destination.unlink()
                elif destination.exists():
                    raise FileExistsError(destination)
                destination.symlink_to(origin)
        combined.append(summary)
        inventory.append(
            {
                "sample_id": sample,
                "species_current": row["species_current"],
                "study": row["study"],
                "source_root": row["result_root"],
                "selection_reason": row["selection_reason"],
                "assembly_status": summary["assembly_status"],
                "assembly_fasta": str(fasta),
                "assembly_sha256": checksum,
                "assembly_bases": summary["assembly_bases"],
                "raw_bins": summary["raw_bins"],
                "viral_contigs": summary["viral_contigs"],
            }
        )
        mapping.append(mapping_input(row, settings))
    write_table(
        root / "phase3_cohort/assembly_summary.tsv", combined, list(combined[0])
    )
    write_table(
        root / "phase3_cohort/membership.tsv",
        cohort,
        [
            "sample_id",
            "host_route",
            "strategy",
            "retained_pairs",
            "minimum_assembly_input_pairs",
            "assembly_eligible",
        ],
    )
    write_table(root / "assembly_sources.tsv", inventory, list(inventory[0]))
    write_table(root / "mapping_inputs.tsv", mapping, list(mapping[0]))
    (root / "assembly_preparation.json").write_text(
        json.dumps(
            {
                "status": "staged_for_qc",
                "mapping_authorized": False,
                "samples": len(cohort),
                "assembly_statuses": dict(
                    Counter(r["assembly_status"] for r in combined)
                ),
                "mapping_input_availability": dict(
                    Counter(r["availability"] for r in mapping)
                ),
                "source_summary_sha256": {
                    str(Path(s) / "phase3_cohort/assembly_summary.tsv"): sha256(
                        Path(s) / "phase3_cohort/assembly_summary.tsv"
                    )
                    for s in summaries
                },
            },
            indent=2,
        )
        + "\n"
    )


def fasta_metrics(path: Path) -> list[dict]:
    records = []
    name, sequence = None, []

    def finish() -> None:
        if name is None:
            return
        seq = "".join(sequence).upper()
        if not seq or set(seq) - set("ACGTN"):
            raise ValueError(f"Invalid contig sequence: {path} {name}")
        acgt = len(seq) - seq.count("N")
        records.append(
            {
                "contig_id": name,
                "length": len(seq),
                "gc_fraction": (seq.count("G") + seq.count("C")) / acgt
                if acgt
                else "NA",
                "n_fraction": seq.count("N") / len(seq),
            }
        )

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            finish()
            name, sequence = line[1:].split()[0], []
        else:
            if name is None and line.strip():
                raise ValueError(f"Sequence before FASTA header: {path}")
            sequence.append(line.strip())
    finish()
    if not records or len({r["contig_id"] for r in records}) != len(records):
        raise ValueError(f"Empty or duplicate-ID candidate FASTA: {path}")
    return records


def review(settings: dict) -> None:
    check_hold()
    root = Path(settings["preparation_root"])
    cohort = unique(table(settings["cohort"]), "sample_id")
    candidates = table(root / "qc/candidates.tsv")
    unique(candidates, "mag_id")
    known = {
        r["bin_id"]: r["sha256"]
        for r in json.loads(Path("config/mag_domain_audit.json").read_text())["bins"]
    }
    excluded = {
        b
        for r in table("config/catalog_domain_exclusions.tsv")
        for b in r["bin_ids"].split(";")
    }
    if set(known) != excluded:
        raise ValueError("Known domain exclusions must have sequence identities")
    old = unique(table("data/results/phase3_catalog/mags/mag_catalog.tsv"), "mag_id")
    rows, contigs, seen = [], [], set()
    for candidate in candidates:
        name, sample = candidate["mag_id"], candidate["sample_id"]
        path = root / "qc/candidates" / f"{name}.fa"
        digest = sha256(path)
        if name in known and digest != known[name]:
            raise ValueError(f"Known eukaryotic control changed sequence: {name}")
        records = fasta_metrics(path)
        if seen & {r["contig_id"] for r in records}:
            raise ValueError(f"Contig appears in multiple candidate bins: {name}")
        seen.update(r["contig_id"] for r in records)
        contigs.extend(
            {"mag_id": name, "sample_id": sample, **r, "gene_domain_qc": "pending"}
            for r in records
        )
        prior = Path("data/results/phase3_catalog/mags/candidates") / f"{name}.fa"
        matched = prior.is_file() and sha256(prior) == digest
        row = {
            **candidate,
            "species_current": cohort[sample]["species_current"],
            "source_root": cohort[sample]["result_root"],
            "fasta": str(path),
            "sha256": digest,
            "contigs": len(records),
            "previous_gtdb_classification": old.get(name, {}).get(
                "gtdb_classification", ""
            )
            if matched
            else "",
            "previous_sequence_matches": str(matched).lower(),
            "known_domain": "Eukaryota" if name in excluded else "",
            "catalog_status": "excluded_known_eukaryote_control"
            if name in excluded
            else "pending_domain_and_contamination_qc",
            "gene_domain_qc": "pending",
            "contig_taxonomy_qc": "pending",
            "review_decision": "pending",
        }
        rows.append(row)
    if not excluded <= {r["mag_id"] for r in rows}:
        raise ValueError(
            "A known eukaryotic positive control is absent from the QC candidates"
        )
    out = Path("figures/catalog_qc")
    write_table(out / "candidate_review.tsv", rows, list(rows[0]))
    write_table(root / "qc/contig_review.tsv", contigs, list(contigs[0]))
    staged = json.loads((root / "assembly_preparation.json").read_text())
    report = {
        **staged,
        "generated_utc": datetime.now(UTC).isoformat(),
        "candidate_bins": len(rows),
        "known_eukaryotic_controls": len(excluded),
        "candidates_requiring_domain_qc": len(rows) - len(excluded),
        "catalog_approved": False,
        "mapping_authorized": False,
        "pending_reviews": [
            "Broad-domain gene classification of every candidate, including existing bins",
            "Contig-level taxonomic coherence and host/eukaryotic contamination",
            "GTDB marker support and classification conflicts",
            "CheckM2 quality in light of verified domain identity",
            "Viral hallmarks, CheckV, and endogenous/host-sequence evidence",
            "Final deduplication and catalog review before mapping authorization",
        ],
        "input_sha256": {
            p: sha256(p)
            for p in [
                settings["cohort"],
                "config/assembly_completion.json",
                "config/phase3_catalog.json",
                "config/catalog_domain_exclusions.tsv",
                "config/post_expansion_remapping.json",
                "scripts/prepare_mapping_review.py",
                str(root / "qc/candidates.tsv"),
            ]
        },
    }
    (out / "preparation_summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: report[k]
                for k in (
                    "samples",
                    "candidate_bins",
                    "known_eukaryotic_controls",
                    "mapping_authorized",
                    "catalog_approved",
                )
            },
            indent=2,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["stage", "review"])
    parser.add_argument(
        "--config", default="config/assembly_completion.json", type=Path
    )
    args = parser.parse_args()
    settings = json.loads(args.config.read_text())
    (stage if args.stage == "stage" else review)(settings)


if __name__ == "__main__":
    main()
