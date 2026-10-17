#!/usr/bin/env python3
"""Export curated publication metadata without changing frozen analysis inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INPUTS = [
    "manifest.csv",
    "config/samples.tsv",
    "config/publication_host_names.tsv",
    "data/metadata/sample_metadata_updates.tsv",
    "data/metadata/publication_metadata_updates.tsv",
]
PUBLIC_FIELDS = [
    "library_id",
    "analysis_sample_id",
    "specimen_id",
    "voucher_id",
    "biosample",
    "collection_id",
    "sample_id",
    "tissue_tube_id",
    "study",
    "study_memberships",
    "species_current",
    "species_as_published",
    "ocean_region",
    "location",
    "collection_site",
    "locality",
    "latitude",
    "longitude",
    "collection_date",
    "depth_m",
    "depth_min_m",
    "depth_max_m",
    "depth_original",
    "collection_depth_source",
    "collection_depth_basis",
    "collection_method",
    "tissue",
    "tissue_status",
    "life_stage",
    "country",
    "specimen_voucher",
    "sra_run",
    "bioproject",
    "sequencing_batches",
    "read_pairs",
    "include_primary",
    "analysis_notes",
    "notes",
]
EDITABLE = {
    "voucher_id",
    "biosample",
    "collection_id",
    "sample_id",
    "tissue_tube_id",
    "location",
    "collection_site",
    "locality",
    "latitude",
    "longitude",
    "collection_date",
    "depth_m",
    "depth_min_m",
    "depth_max_m",
    "depth_original",
    "collection_depth_source",
    "collection_depth_basis",
    "collection_method",
    "tissue",
    "life_stage",
    "country",
    "sra_run",
    "notes",
    "analysis_notes",
}
AUDIT_FIELDS = ["library_id", "field", "value", "source", "source_locator", "notes"]
CHECK_FIELDS = [
    "voucher_id",
    "biosample",
    "collection_id",
    "sample_id",
    "latitude",
    "longitude",
    "depth_m",
    "location",
    "tissue",
]
# Provenance also covers optional fields added during museum-record review.
# The missing-field checklist remains focused on the required sample metadata.
PROVENANCE_FIELDS = [
    *CHECK_FIELDS,
    "collection_site",
    "collection_date",
    "collection_method",
    "country",
    "depth_min_m",
    "depth_max_m",
    "depth_original",
    "collection_depth_source",
    "collection_depth_basis",
    "tissue_tube_id",
    "life_stage",
    "notes",
]
MISSING = {
    "",
    "na",
    "n/a",
    "missing",
    "not provided",
    "not collected",
    "not applicable",
}


def missing(value: str) -> bool:
    """Recognize source missing-value labels without treating numeric zero as missing."""
    return value.strip().lower() in MISSING or value.startswith("TODO:")


def human_notes(row: dict[str, str]) -> str:
    """Explain missing evidence and derived values in specimen-facing language."""
    notes = [row["notes"]] if row.get("notes") else []
    descriptions = {
        "voucher_id": "No museum voucher number was found in the reviewed source records.",
        "biosample": "No BioSample accession was identified in the reviewed study and NCBI records.",
        "collection_id": "Field collection-event identifier not recorded in the reviewed source tables.",
        "sample_id": "No separate tissue/sample label was found in the reviewed source tables.",
    }
    for field, description in descriptions.items():
        if missing(row.get(field, "")):
            notes.append(description)
    if missing(row.get("tissue", "")):
        notes.append(
            "Tissue is unspecified in the NCBI sample record and the reviewed study tables."
            if not missing(row.get("biosample", ""))
            else "Tissue was not specified in the reviewed source records."
        )
    basis = row.get("collection_depth_basis", "")
    if basis == "interval_midpoint":
        notes.append(
            f"Depth is the midpoint of the reported {row['depth_min_m']}–{row['depth_max_m']} m range."
        )
    elif basis == "curator_assigned_surface":
        notes.append(
            "Depth assigned as 0 m under the Physalia surface convention; no specimen-specific measured depth was recorded."
        )
    if row.get("tissue_tube_id") and row.get("sample_id") == row["tissue_tube_id"]:
        notes.append("Sample ID is the tissue-tube label recorded by the source study.")
    return " ".join(notes)


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(
            csv.DictReader(handle, delimiter="," if path.suffix == ".csv" else "\t")
        )


def unique_index(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    index = {}
    for row in rows:
        value = row[key]
        if not value or value in index:
            raise ValueError(f"Empty or duplicate {key}: {value!r}")
        index[value] = row
    return index


def curate_metadata(
    manifest: list[dict[str, str]],
    samples: list[dict[str, str]],
    host_names: list[dict[str, str]],
    updates: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Join immutable workflow handles to biological identifiers and reviewed values."""
    original = unique_index(manifest, "library_id")
    workflow = unique_index(samples, "library_id")
    unique_index(samples, "sample_id")
    if set(original) != set(workflow):
        raise ValueError("Manifest and workflow library membership differ")
    rows = {}
    for library, source in original.items():
        if source["specimen_id"] != workflow[library]["specimen_id"]:
            raise ValueError(f"Workflow specimen_id differs: {library}")
        rows[library] = {
            **source,
            "analysis_sample_id": workflow[library]["sample_id"],
            "voucher_id": source.get("specimen_voucher", ""),
            "sample_id": "",
            "location": source.get("locality", ""),
            "collection_site": "",
            "analysis_notes": source.get("notes", ""),
            "notes": "",
        }
    baseline = {library: dict(row) for library, row in rows.items()}
    applied = {}
    for update in updates:
        library, field, value = (update[k] for k in ("library_id", "field", "value"))
        if library not in rows or field not in EDITABLE:
            raise ValueError(f"Unknown library or protected field: {library}, {field}")
        if not update["source"] or not update["source_locator"]:
            raise ValueError(f"Missing provenance: {library}, {field}")
        if update["original_value"] != baseline[library].get(field, ""):
            raise ValueError(f"Stale original_value: {library}, {field}")
        key = (library, field)
        if key in applied and applied[key] != value:
            raise ValueError(f"Conflicting publication updates: {library}, {field}")
        applied[key] = value
        rows[library][field] = value

    for library, update in unique_index(host_names, "library_id").items():
        if library not in rows:
            raise ValueError(f"Unknown host-name update: {library}")
        row = rows[library]
        if row["species_current"] != update["previous_name"]:
            raise ValueError(f"Unexpected previous host name: {library}")
        if update["assignment"] not in row["analysis_notes"].split("; "):
            raise ValueError(f"Host assignment not supported: {library}")
        row["species_current"] = update["publication_name"]

    for library, row in rows.items():
        for field, limit in (("latitude", 90), ("longitude", 180)):
            if not missing(row.get(field, "")):
                value = Decimal(row[field])
                if not value.is_finite() or abs(value) > limit:
                    raise ValueError(f"Invalid {field}: {library}")
        depths = {}
        for field in ("depth_m", "depth_min_m", "depth_max_m"):
            if not missing(row.get(field, "")):
                value = Decimal(row[field])
                if not value.is_finite() or value < 0:
                    raise ValueError(f"Invalid {field}: {library}")
                depths[field] = value
        if "depth_min_m" in depths and "depth_max_m" in depths:
            low, high = depths["depth_min_m"], depths["depth_max_m"]
            if low > high:
                raise ValueError(f"Reversed depth interval: {library}")
            if "depth_m" in depths and not low <= depths["depth_m"] <= high:
                raise ValueError(f"Point depth outside interval: {library}")
            if "depth_m" not in depths:
                row["depth_m"] = format((low + high) / 2, "f")
                row["collection_depth_basis"] = "interval_midpoint"
        row["tissue_status"] = (
            "not_recorded" if missing(row.get("tissue", "")) else "recorded"
        )
        # Retain the older column as an exact alias for existing table consumers.
        row["specimen_voucher"] = row["voucher_id"]
        row["notes"] = human_notes(row)
    return list(rows.values())


def load_metadata(root: Path = ROOT) -> list[dict[str, str]]:
    """Shared entry point for the metadata-only export and cohort figure builder."""
    return curate_metadata(
        read_table(root / INPUTS[0]),
        read_table(root / INPUTS[1]),
        read_table(root / INPUTS[2]),
        read_table(root / INPUTS[4]),
    )


def write_table(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    """Avoid touching an unchanged export, including its modification time."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    writer.writerows(rows)
    content = buffer.getvalue().encode()
    if not path.exists() or path.read_bytes() != content:
        path.write_bytes(content)


def export_metadata(output_dir: Path) -> list[dict[str, str]]:
    """Write metadata and its audit together, also when rebuilding cohort figures."""
    rows = load_metadata()
    sources = defaultdict(list)
    for path in (INPUTS[3], INPUTS[4]):
        for record in read_table(ROOT / path):
            sources[(record["library_id"], record["field"])].append(record)
    audit, gaps = [], []
    for row in rows:
        for field in PROVENANCE_FIELDS:
            value = row.get(field, "")
            evidence = sources.get((row["library_id"], field), [])
            if field == "notes":
                evidence = [
                    *evidence,
                    {
                        "source": "scripts/export_sample_metadata.py",
                        "source_locator": f"human_notes;library_id={row['library_id']}",
                        "notes": "Generated explanations of missing fields and depth conventions supplement any sourced notes.",
                    },
                ]
            if (
                field == "depth_m"
                and row.get("collection_depth_basis") == "interval_midpoint"
            ):
                evidence = [
                    {
                        "source": row.get("collection_depth_source") or "manifest.csv",
                        "source_locator": "depth_min_m;depth_max_m",
                        "notes": "User convention: (minimum + maximum) / 2; interval preserved.",
                    }
                ]
            if not evidence:
                evidence = [
                    {
                        "source": "manifest.csv",
                        "source_locator": f"library_id={row['library_id']};{field}",
                        "notes": "Frozen baseline; empty values are unresolved, not inferred.",
                    }
                ]
            for record in evidence:
                audit.append(
                    {
                        "library_id": row["library_id"],
                        "field": field,
                        "value": value,
                        **{
                            k: record.get(k, "")
                            for k in ("source", "source_locator", "notes")
                        },
                    }
                )
            if field in CHECK_FIELDS and missing(value):
                gaps.append(
                    {
                        "library_id": row["library_id"],
                        "analysis_sample_id": row["analysis_sample_id"],
                        "field": field,
                        "status": "not_recorded_or_unresolved",
                    }
                )
    output_dir.mkdir(parents=True, exist_ok=True)
    write_table(output_dir / "sample_metadata.tsv", rows, PUBLIC_FIELDS)
    write_table(output_dir / "sample_metadata_provenance.tsv", audit, AUDIT_FIELDS)
    write_table(
        output_dir / "sample_metadata_missing.tsv",
        gaps,
        ["library_id", "analysis_sample_id", "field", "status"],
    )
    record = {
        "command": [sys.executable, *sys.argv],
        "python": sys.version.split()[0],
        "libraries": len(rows),
        "missing": dict(Counter(row["field"] for row in gaps)),
        "depth_basis": dict(
            Counter(row.get("collection_depth_basis", "") for row in rows)
        ),
        "inputs_sha256": {
            p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()
            for p in [*INPUTS, "scripts/export_sample_metadata.py"]
        },
        "outputs_sha256": {
            name: hashlib.sha256((output_dir / name).read_bytes()).hexdigest()
            for name in [
                "sample_metadata.tsv",
                "sample_metadata_provenance.tsv",
                "sample_metadata_missing.tsv",
            ]
        },
    }
    (output_dir / "sample_metadata_export.json").write_text(
        json.dumps(record, indent=2) + "\n"
    )
    print(
        json.dumps(
            {key: record[key] for key in ("libraries", "missing", "depth_basis")},
            indent=2,
        )
    )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/cohort")
    args = parser.parse_args()
    export_metadata(args.output_dir)


if __name__ == "__main__":
    main()
