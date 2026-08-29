#!/usr/bin/env python3
"""Validate the frozen Phase-0 publication input tables.

This performs metadata, identity, routing, and filesystem checks without
streaming every full FASTQ. Full read counts come from SRA spots or the recorded
completed count jobs; FASTQ parsing and pair synchronization are rechecked by
the preprocessing workflow.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import os
from collections import Counter, defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = REPO / "manifest.csv"
DEFAULT_RAW = REPO / "data" / "metadata" / "raw_files.tsv"
DEFAULT_PROV = REPO / "data" / "metadata" / "library_provenance.tsv"
DEFAULT_RESOURCES = REPO / "data" / "metadata" / "input_resources.tsv"
DEFAULT_AUDITS = REPO / "data" / "metadata" / "reference_audits.tsv"
DEFAULT_FREEZE = REPO / "data" / "metadata" / "manifest.freeze.sha256"
DEFAULT_REPORT = REPO / "data" / "metadata" / "manifest_validation.txt"


def load(path: Path, delimiter: str) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def split(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def validate(manifest_path: Path, raw_path: Path, prov_path: Path,
             resources_path: Path, audits_path: Path,
             freeze_path: Path | None) -> list[str]:
    errors: list[str] = []

    def require(condition: bool, message: str) -> None:
        if not condition:
            errors.append(message)

    manifest = load(manifest_path, ",")
    raw = load(raw_path, "\t")
    provenance = load(prov_path, "\t")
    resources = load(resources_path, "\t")
    audits = load(audits_path, "\t")
    ids = [row["library_id"] for row in manifest]
    specimens = [row["specimen_id"] for row in manifest]
    by_id = {row["library_id"]: row for row in manifest}

    require(len(manifest) == 205, f"expected 205 manifest rows, found {len(manifest)}")
    require(len(set(ids)) == len(ids), "library_id values are not unique")
    require(len(set(specimens)) == len(specimens), "specimen_id values are not unique")
    require(Counter(row["study"] for row in manifest) == Counter({
        "Church2025": 151, "Ahuja2024": 33, "Ahuja2026": 21,
    }), "unexpected primary-study counts")
    require(all(row["include_primary"] == "true" for row in manifest),
            "not all manifest rows are included in the primary analysis")
    require(all(not row["exclusion_reason"] for row in manifest),
            "included rows have exclusion reasons")
    require(not any(
        value.startswith("TODO") for row in manifest for value in row.values()
    ), "manifest contains unresolved TODO values")
    require(all(int(row["read_pairs"]) > 0 for row in manifest),
            "read_pairs must be positive integers")
    require(all(row["depth_source"] in {"SRA_spots", "counted"} for row in manifest),
            "unexpected depth_source")
    require(all(row["host_reference"] in {"P_physalis", "N_septata", "none"}
                for row in manifest), "unexpected host_reference")
    for row in manifest:
        species = row["species_current"]
        if species.startswith("Physalia"):
            require(row["host_reference"] == "P_physalis",
                    f"Physalia routing error: {row['library_id']}")
        if species == "Nanomia septata":
            require(row["host_reference"] == "N_septata",
                    f"N. septata routing error: {row['library_id']}")

    raw_by_library: dict[str, list[dict[str, str]]] = defaultdict(list)
    raw_keys: list[tuple[str, str]] = []
    require(len(raw) == 313, f"expected 313 paired FASTQ rows, found {len(raw)}")
    for row in raw:
        library_id = row["library_id"]
        require(library_id in by_id, f"raw row has unknown library: {library_id}")
        raw_by_library[library_id].append(row)
        r1, r2 = Path(row["r1_path"]), Path(row["r2_path"])
        require(r1.is_file(), f"missing R1: {r1}")
        require(r2.is_file(), f"missing R2: {r2}")
        if r1.is_file():
            require(r1.stat().st_size == int(row["r1_bytes"]),
                    f"R1 size changed: {r1}")
            with r1.open("rb") as fh:
                require(fh.read(2) == b"\x1f\x8b", f"R1 lacks gzip magic: {r1}")
        if r2.is_file():
            require(r2.stat().st_size == int(row["r2_bytes"]),
                    f"R2 size changed: {r2}")
            with r2.open("rb") as fh:
                require(fh.read(2) == b"\x1f\x8b", f"R2 lacks gzip magic: {r2}")
        require(int(row["r1_read_length"]) > 0 and int(row["r2_read_length"]) > 0,
                f"invalid read length: {library_id}/{row['read_pair_id']}")
        require("combined_R" not in r1.name and "combined_R" not in r2.name,
                f"derived combined FASTQ selected: {library_id}")
        require(row["sequencing_batch"] and "lane_unknown" not in row["sequencing_batch"],
                f"unresolved sequencing batch: {library_id}/{row['read_pair_id']}")
        require(row["sequencing_batch_source"] ==
                "FASTQ header: instrument/run/flowcell/lane",
                f"unexpected batch source: {library_id}/{row['read_pair_id']}")
        raw_keys.append((os.path.realpath(r1), os.path.realpath(r2)))
    require(len(set(raw_keys)) == len(raw_keys),
            "a paired FASTQ set is assigned to more than one analytical library")
    require(set(raw_by_library) == set(ids),
            "manifest and raw-files tables have different library sets")

    for library_id, row in by_id.items():
        rows = raw_by_library[library_id]
        require(len(rows) == int(row["n_lanes"]),
                f"n_lanes mismatch: {library_id}")
        require(set(split(row["r1_paths"])) == {value["r1_path"] for value in rows},
                f"R1 aggregate mismatch: {library_id}")
        require(set(split(row["r2_paths"])) == {value["r2_path"] for value in rows},
                f"R2 aggregate mismatch: {library_id}")
        require(set(split(row["sequencing_batches"])) ==
                {value["sequencing_batch"] for value in rows},
                f"batch aggregate mismatch: {library_id}")
        require(set(split(row["raw_path_mccleary"])) ==
                {value["source_directory"] for value in rows},
                f"source-directory aggregate mismatch: {library_id}")

    require(len(provenance) == 208,
            f"expected 208 provenance rows, found {len(provenance)}")
    prov_ids = [row["provenance_id"] for row in provenance]
    require(len(set(prov_ids)) == len(prov_ids), "provenance_id values are not unique")
    require(all(row["library_id"] in by_id for row in provenance),
            "provenance references unknown analytical library")
    prov_by_library: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in provenance:
        prov_by_library[row["library_id"]].append(row)
        if row["library_id"] in by_id:
            require(row["specimen_id"] == by_id[row["library_id"]]["specimen_id"],
                    f"provenance specimen mismatch: {row['provenance_id']}")
    require(set(prov_by_library) == set(ids),
            "manifest and provenance tables have different library sets")
    for library_id, row in by_id.items():
        studies = {value["study"] for value in prov_by_library[library_id]}
        require(set(split(row["study_memberships"])) == studies,
                f"study-membership mismatch: {library_id}")
        require(set(split(row["also_in_studies"])) == studies - {row["study"]},
                f"also_in_studies mismatch: {library_id}")

    # Explicit identity/provenance invariants behind the Phase-0 reconciliation.
    require("Ahuja2024:NA22" not in by_id, "NA22 duplicate remains analytical")
    require(by_id.get("Church2025:YPM-IZ-104465", {}).get("study_memberships") ==
            "Church2025;Ahuja2024", "NA22/Church provenance was not merged")
    require(by_id.get("Ahuja2024:CWD16", {}).get("specimen_id") == "YPM:IZ:35039",
            "CWD16 identity changed")
    require(by_id.get("Ahuja2026:NA19", {}).get("specimen_id") == "YPM:IZ:111744",
            "NA19 identity changed")
    require(by_id.get("Ahuja2024:CWD16", {}).get("study_memberships") ==
            "Ahuja2024;Ahuja2026", "CWD16 reuse provenance missing")
    require(by_id.get("Ahuja2026:NA19", {}).get("study_memberships") ==
            "Ahuja2024;Ahuja2026", "NA19 origin provenance missing")

    require(len(resources) == 10,
            f"expected 10 Phase-0 resources, found {len(resources)}")
    require(len({row["resource_id"] for row in resources}) == len(resources),
            "resource_id values are not unique")
    for row in resources:
        path, anchor = Path(row["path"]), Path(row["checksum_anchor"])
        require(row["status"] == "ready", f"resource not ready: {row['resource_id']}")
        require(path.exists(), f"resource path missing: {path}")
        require(anchor.is_file(), f"resource checksum anchor missing: {anchor}")
        if anchor.is_file():
            require(sha256(anchor) == row["anchor_sha256"],
                    f"resource checksum changed: {row['resource_id']}")

    require(len(audits) == 2, f"expected 2 host audits, found {len(audits)}")
    require({row["reference_accession"] for row in audits} ==
            {"GCA_041430235.2", "GCA_048301705.1"},
            "host-reference audit accessions are incomplete")
    for row in audits:
        root = Path(row["evidence_root"])
        fai = root.parent / "unplaced.fasta.fai"
        genes = root / "unplaced_annotate" / "unplaced_genes.tsv"
        proviruses = root / "unplaced_find_proviruses" / "unplaced_provirus.tsv"
        require(root.is_dir(), f"host-audit evidence root missing: {root}")
        require(fai.is_file(), f"host-audit FAI missing: {fai}")
        require(genes.is_file(), f"host-audit genes missing: {genes}")
        require(proviruses.is_file(), f"host-audit provirus table missing: {proviruses}")
        if fai.is_file():
            entries = [line.split("\t") for line in fai.read_text().splitlines()]
            require(len(entries) == int(row["unplaced_scaffolds"]),
                    f"host-audit scaffold count changed: {row['reference_accession']}")
            require(sum(int(entry[1]) for entry in entries) == int(row["unplaced_bases"]),
                    f"host-audit base count changed: {row['reference_accession']}")
        if genes.is_file():
            gene_rows = load(genes, "\t")
            require(sum(int(value["plasmid_hallmark"]) for value in gene_rows) ==
                    int(row["plasmid_hallmark_genes"]),
                    f"plasmid hallmark count changed: {row['reference_accession']}")
            require(sum(value["annotation_conjscan"] != "NA" for value in gene_rows) ==
                    int(row["conjugation_genes"]),
                    f"conjugation count changed: {row['reference_accession']}")
        if proviruses.is_file():
            require(len(load(proviruses, "\t")) == int(row["provirus_calls"]),
                    f"provirus count changed: {row['reference_accession']}")

    if freeze_path is not None:
        require(freeze_path.is_file(), f"missing freeze file: {freeze_path}")
        if freeze_path.is_file():
            for line in freeze_path.read_text().splitlines():
                expected, relative = line.split(None, 1)
                path = REPO / relative.strip()
                require(path.is_file(), f"frozen input missing: {relative.strip()}")
                if path.is_file():
                    require(sha256(path) == expected,
                            f"checksum mismatch: {relative.strip()}")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--raw-files", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--provenance", type=Path, default=DEFAULT_PROV)
    parser.add_argument("--resources", type=Path, default=DEFAULT_RESOURCES)
    parser.add_argument("--reference-audits", type=Path, default=DEFAULT_AUDITS)
    parser.add_argument("--freeze", type=Path, default=None,
                        help="verify a freeze checksum file")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    errors = validate(args.manifest, args.raw_files, args.provenance,
                      args.resources, args.reference_audits, args.freeze)
    status = "PASS" if not errors else "FAIL"
    lines = [
        f"status\t{status}",
        "manifest_libraries\t205",
        "unique_specimens\t205",
        "paired_fastq_records\t313",
        "provenance_records\t208",
        "phase0_resources\t10",
        "host_reference_audits\t2",
        f"freeze_checked\t{'yes' if args.freeze else 'no'}",
        f"errors\t{len(errors)}",
    ] + [f"error\t{error}" for error in errors]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
