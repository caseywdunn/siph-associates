#!/usr/bin/env python3
"""Package curated sequences with source metadata and explicit submission blockers."""

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from prepare_deposit_mag import read_fasta, read_table


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset:offset + 80] + "\n")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    base = Path(config["output_root"])
    bundle = base / "bundle"
    decisions = read_table(Path(config["decisions"]))
    names = {r["mag_id"]: r for r in read_table(Path(config["organism_names"]))}
    catalog = {r["mag_id"]: r for r in read_table(Path(config["mag_catalog"])) if r["catalog_status"] == "retained"}
    with Path(config["sample_manifest"]).open(newline="") as handle:
        samples = {r["library_id"].replace(":", "__"): r for r in csv.DictReader(handle)}
    manifest, metadata, sequence_map, quality = [], [], [], []
    sets = {"mags": [], "viruses": [], "viral_endogenous": [], "ssu": [], "markers": []}
    viral_manifests = {
        "viruses": read_table(Path(config["other_sequences"]["viruses"]["source_table"])),
        "viral_endogenous": read_table(Path(config["viral_endogenous_manifest"])),
    }
    viral_rows = {}
    viral_classes = {}
    for group, rows in viral_manifests.items():
        for row in rows:
            identifier = row["votu_id"]
            if identifier in viral_rows:
                raise ValueError(f"Repeated viral catalog identifier across reporting classes: {identifier}")
            viral_rows[identifier] = row
            viral_classes[identifier] = group
    for row in decisions:
        mag = row["mag_id"]
        genome = catalog[mag]
        sample = samples[genome["sample_id"]]
        taxonomy = names[mag]
        edit = json.loads((base / "edits" / f"{mag}.json").read_text())
        if edit["edits"]:
            report_path = base / "checkm2" / mag / "quality_report.tsv"
            report = read_table(report_path)
            if len(report) != 1 or report[0]["Name"] != mag:
                raise ValueError(f"CheckM2 report does not match edited MAG: {mag}")
            complete, contamination = float(report[0]["Completeness"]), float(report[0]["Contamination"])
            quality_source = str(report_path)
        else:
            complete, contamination = float(genome["completeness"]), float(genome["contamination"])
            quality_source = config["mag_catalog"]
        passes_numeric_thresholds = complete >= 90 and edit["genome_size"] >= 100000
        exclusions = []
        holds = []
        if complete < 90:
            exclusions.append("completeness_below_90_percent")
        if edit["genome_size"] < 100000:
            exclusions.append("genome_below_100kb")
        if row["deposit_exclusion"]:
            holds.append(row["deposit_exclusion"])
        if edit["remaining_exact_palindromes"]:
            holds.append("unreviewed_exact_palindrome")
        route = "held_provenance_only" if holds else "archive" if exclusions else "genbank_candidate"
        exclusions.extend(holds)
        accessions = [a for a in re.split(r"[;, ]+", sample["sra_run"]) if a]
        if any(not re.fullmatch(r"[SED]RR\d+", a) for a in accessions):
            raise ValueError(f"Invalid source-read accession: {mag}: {accessions}")
        blockers = list(exclusions)
        if not accessions:
            blockers.append("source_read_accessions_pending")
        if taxonomy["name_status"] != "verified":
            blockers.append("ncbi_organism_name_verification_or_creation_pending")
        blockers += ["author_ncbi_name_confirmation_pending", "mimag_biosample_registration_pending",
                     "fcs_gx_pending", "table2asn_pending", "pgap_taxcheck_pending"]
        records = []
        for index, (original, sequence) in enumerate(read_fasta(base / "cleaned_mags" / f"{mag}.fa").items(), 1):
            sequence_id = f"{taxonomy['isolate']}_{index:05d}"
            header = sequence_id + (f" [SRA={','.join(accessions)}]" if accessions else "")
            records.append((header, sequence))
            sequence_map.append({"set": "mags", "mag_id": mag, "original_id": original,
                                 "submission_id": sequence_id, "length": len(sequence)})
        path = bundle / "mags" / f"{mag}.fna"
        digest = write_fasta(path, records)
        sets["mags"].extend(records)
        manifest.append({"mag_id": mag, "species_cluster": row["species_cluster"], "representative": row["representative"],
                         "route": route, "passes_numeric_thresholds": str(passes_numeric_thresholds).lower(),
                         "sequence_review_status": "held_pending_review" if holds else "reviewed_edits_applied", "fasta": str(path), "sha256": digest, "genome_size": edit["genome_size"],
                         "contigs": edit["contigs"], "completeness": complete, "contamination": contamination,
                         "quality_source": quality_source, "proposed_ncbi_name": taxonomy["proposed_ncbi_name"],
                         "ncbi_taxid": taxonomy["ncbi_taxid"], "isolate": taxonomy["isolate"],
                         "source_reads": ",".join(accessions), "submission_status": "not_submitted",
                         "blockers": ";".join(blockers)})
        metadata.append({"sample_name": taxonomy["isolate"], "package": "MIMAG Metagenome-assembled Genome",
                         "organism": taxonomy["proposed_ncbi_name"], "taxid": taxonomy["ncbi_taxid"],
                         "isolate": taxonomy["isolate"], "sample_type": "metagenomic assembly",
                         "derived_from": ",".join(accessions), "host": sample["species_current"],
                         "geo_loc_name": sample["country"], "collection_date": sample["collection_date"],
                         "lat_lon": sample["lat_long_raw"], "isolation_source": sample["tissue"],
                         "source_biosample": sample["biosample"], "MAG_biosample": "", "MAG_bioproject": "",
                         "completeness_score": complete, "contamination_score": contamination,
                         "completion_check": "CheckM2", "metadata_status": "draft_author_confirmation_required"})
        updated = dict(genome)
        updated.update(completeness=complete, contamination=contamination,
                       genome_size=edit["genome_size"], gc_content=edit["gc_fraction"],
                       quality_score=complete - 5 * contamination,
                       near_complete=str(complete >= 90 and contamination < 5).lower(),
                       sequence_version="curated_for_deposit", quality_source=quality_source,
                       original_catalog_sequence_unchanged="true")
        quality.append(updated)
    for label, spec in config["other_sequences"].items():
        sequences = read_fasta(Path(spec["fasta"]))
        if label == "viruses" and set(sequences) != set(viral_rows):
            raise ValueError("Full viral FASTA IDs do not match the union of the current reporting manifests")
        for index, (original, sequence) in enumerate(sequences.items(), 1):
            if label == "viruses" and hashlib.sha256(sequence.encode()).hexdigest() != viral_rows[original]["sequence_sha256"]:
                raise ValueError(f"Viral FASTA sequence hash differs from reporting manifest: {original}")
            sequence_id = f"SIPH_{label}_{index:05d}"
            destination_label = viral_classes[original] if label == "viruses" else label
            sets[destination_label].append((sequence_id, sequence))
            sequence_map.append({"set": destination_label, "mag_id": "", "original_id": original,
                                 "submission_id": sequence_id, "length": len(sequence)})
        # Preserve source and feature evidence without inventing organism names or CDS boundaries.
        source_path = Path(spec["source_table"])
        destination = bundle / "source_evidence" / source_path.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(source_path.read_bytes())
    source_path = Path(config["viral_endogenous_manifest"])
    (bundle / "source_evidence" / source_path.name).write_bytes(source_path.read_bytes())
    if any(len(sets[group]) != len(rows) for group, rows in viral_manifests.items()):
        raise ValueError("Packaged viral class counts differ from the current reporting manifests")
    set_manifest = []
    for label, records in sets.items():
        path = bundle / "sets" / f"{label}.fna"
        digest = write_fasta(path, records)
        set_manifest.append({"set": label, "sequences": len(records), "fasta": str(path), "sha256": digest,
                             "fcs_adaptor_status": "pending", "feature_validation_status": "pending"})
    write_table(bundle / "mag_manifest.tsv", manifest)
    write_table(bundle / "mimag_biosample_metadata.tsv", metadata)
    write_table(bundle / "sequence_id_mapping.tsv", sequence_map)
    write_table(bundle / "sequence_sets.tsv", set_manifest)
    write_table(bundle / "curated_mag_catalog.tsv", quality)
    write_table(bundle / "genbank_candidates.tsv", [r for r in manifest if r["route"] == "genbank_candidate"])
    write_table(bundle / "archive_genomes.tsv", [r for r in manifest if r["route"] == "archive"])
    write_table(bundle / "held_genomes.tsv", [r for r in manifest if r["route"] == "held_provenance_only"])
    report = {"MAGs": len(manifest),
              "MAGs_passing_numeric_thresholds": sum(r["passes_numeric_thresholds"] == "true" for r in manifest),
              "genbank_candidates": sum(r["route"] == "genbank_candidate" for r in manifest),
              "archive_genomes": sum(r["route"] == "archive" for r in manifest),
              "held_provenance_only": sum(r["route"] == "held_provenance_only" for r in manifest),
              "next_step": "Run FCS-GX after adaptor review; review reports before final validation.",
              "submission_status": "not_submitted", "sequence_version": "curated; original mapping catalog retained",
              "pending": ["FCS-GX screen and report review", "table2asn", "PGAP taxcheck", "ribovore",
                          "confirm source annotations and COI feature coordinates", "BioSample/BioProject registration"]}
    (bundle / "status.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
