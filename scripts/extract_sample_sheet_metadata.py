#!/usr/bin/env python3
"""Extract auditable metadata candidates without changing the cohort manifest.

Requires openpyxl 3.1.5. Source spreadsheets are preserved in data/sources/
(already ignored as binary source supplements); normalized candidates and cell
provenance are written separately. RNA libraries never supply DNA metadata.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path

import openpyxl

FIELDS = [
    "library_id",
    "field",
    "value",
    "source_type",
    "source",
    "source_locator",
    "source_record",
    "notes",
]
MISSING = {"", "missing", "not collected", "not applicable", "not available"}


def text(value: object) -> str:
    """Preserve spreadsheet values while normalizing dates and whitespace."""
    if value is None:
        return ""
    if isinstance(value, (date, datetime)):
        return (
            value.date().isoformat()
            if isinstance(value, datetime)
            else value.isoformat()
        )
    return str(value).strip()


def voucher(value: str) -> str:
    """Normalize punctuation only, retaining the complete museum identifier."""
    return re.sub(r"[\s:\-]", "", value).upper()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def coordinates(value: str) -> tuple[str, str]:
    """Accept signed decimal pairs or explicit N/S E/W decimal degrees."""
    pair = re.fullmatch(
        r"\s*([+\-]?\d+(?:\.\d+)?)\s*,\s*([+\-]?\d+(?:\.\d+)?)\s*", value
    )
    compass = re.fullmatch(
        r"\s*(\d+(?:\.\d+)?)\s*([NS])\s+(\d+(?:\.\d+)?)\s*([EW])\s*",
        value,
        re.IGNORECASE,
    )
    if pair:
        lat, lon = map(float, pair.groups())
    elif compass:
        a, north, b, east = compass.groups()
        lat = float(a) * (-1 if north.upper() == "S" else 1)
        lon = float(b) * (-1 if east.upper() == "W" else 1)
    else:
        raise ValueError(f"Unrecognized coordinates: {value!r}")
    if not (-90 <= lat <= 90 and -180 <= lon <= 180):
        raise ValueError(f"Coordinates out of bounds: {value!r}")
    return format(lat, ".10g"), format(lon, ".10g")


def depths(value: str) -> tuple[list[tuple[str, str]], str]:
    """Retain intervals; interpret bare numbers only with an explicit note."""
    number = r"\d+(?:\.\d+)?"
    match = re.fullmatch(
        rf"({number})(?:\s*[-–]\s*({number}))?\s*(m|meters?|metres?)?",
        value,
        re.IGNORECASE,
    )
    if not match:
        raise ValueError(f"Unrecognized depth: {value!r}")
    low, high, unit = match.groups()
    note = (
        "explicit_metre_unit"
        if unit
        else "metres_inferred_from_other_depth_cells_in_same_column"
    )
    if high is not None:
        if float(low) > float(high):
            raise ValueError(f"Reversed depth interval: {value}")
        return [("depth_min_m", low), ("depth_max_m", high)], note
    return [("depth_m", low)], note


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path("."))
    args = parser.parse_args()
    repo = args.repo.resolve()
    out = repo / "data/sources/metadata_enrichment"
    out.mkdir(parents=True, exist_ok=True)
    with (repo / "manifest.csv").open() as handle:
        manifest = list(csv.DictReader(handle))
    index: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in manifest:
        library = row["library_id"]
        index[("voucher", voucher(row["specimen_id"]))].add(library)
        for run in row["sra_run"].split(";"):
            if run:
                index[("run", run)].add(library)
    for row in read_tsv(repo / "data/metadata/library_provenance.tsv"):
        index[(row["study"], row["original_label"])].add(row["library_id"])
        if row["sra_run"]:
            index[("run", row["sra_run"])].add(row["library_id"])

    specs = [
        (
            "Supplementary_Table_1.xlsx",
            "Ahuja2024_Supplementary_Table_1.xlsx",
            "Ahuja2024",
            13,
        ),
        (
            "journal.pone.0351247.s008.xlsx",
            "Ahuja2026_journal.pone.0351247.s008.xlsx",
            "Ahuja2026",
            1,
        ),
    ]
    candidates, held, raw, joins, sources = [], [], [], [], []
    for input_name, saved_name, study, header_row in specs:
        source_input = args.input_dir / input_name
        source = repo / "data/sources" / saved_name
        digest = hashlib.sha256(source_input.read_bytes()).hexdigest()
        if (
            source.exists()
            and hashlib.sha256(source.read_bytes()).hexdigest() != digest
        ):
            raise ValueError(
                f"Refusing to replace a different preserved source: {source}"
            )
        if not source.exists():
            shutil.copyfile(source_input, source)
        workbook = openpyxl.load_workbook(source, data_only=True)
        if len(workbook.worksheets) != 1:
            raise ValueError(f"Unexpected additional sheets in {source}")
        sheet = workbook.worksheets[0]
        headers = {
            cell.column: text(cell.value).lstrip("*")
            for cell in sheet[header_row]
            if cell.value
        }
        sources.append(
            {
                "source": str(source.relative_to(repo)),
                "input": str(source_input),
                "sha256": digest,
                "study": study,
                "sheet": sheet.title,
                "header_row": header_row,
                "columns": [
                    {
                        "column": column,
                        "name": name,
                        "comment": sheet.cell(header_row, column).comment.text
                        if sheet.cell(header_row, column).comment
                        else "",
                    }
                    for column, name in headers.items()
                ],
            }
        )
        for number in range(header_row + 1, sheet.max_row + 1):
            values = {
                name: text(sheet.cell(number, column).value)
                for column, name in headers.items()
            }
            if not any(values.values()):
                continue
            cells = {
                name: sheet.cell(number, column).coordinate
                for column, name in headers.items()
            }
            source_record = (
                f"{study}:{values.get('sample_name', values.get('sample_ID', ''))}"
            )
            raw.append(
                {
                    "source": str(source.relative_to(repo)),
                    "sheet": sheet.title,
                    "row": number,
                    "values": values,
                    "cells": cells,
                }
            )
            probes = [
                (study, values.get("sample_name", values.get("sample_ID", ""))),
                (
                    "voucher",
                    voucher(values.get("specimen_voucher", values.get("YPM ID", ""))),
                ),
                ("run", values.get("Accession number", "")),
            ]
            matches = [
                (kind, key, index[(kind, key)])
                for kind, key in probes
                if key and index.get((kind, key))
            ]
            targets = set().union(*(ids for _, _, ids in matches)) if matches else set()
            status = (
                "matched"
                if len(targets) == 1
                else "unmatched"
                if not targets
                else "conflicting_identifiers"
            )
            if study == "Ahuja2026" and values.get("DNA/RNA") != "DNA":
                status = "excluded_non_DNA"
            joined = {
                "source": str(source.relative_to(repo)),
                "sheet": sheet.title,
                "row": str(number),
                "source_record": source_record,
                "status": status,
                "library_id": ";".join(sorted(targets)),
                "matched_identifiers": ";".join(
                    f"{kind}:{key}" for kind, key, _ in matches
                ),
            }
            joins.append(joined)
            if status != "matched":
                continue
            library = next(iter(targets))
            context = {
                "library_id": library,
                "source_type": "primary_sheet",
                "source": str(source.relative_to(repo)),
                "source_record": source_record,
                "notes": joined["matched_identifiers"],
            }
            locators = {key: f"{sheet.title}!{cell}" for key, cell in cells.items()}

            def add(
                field: str,
                value: str,
                original: str,
                note: str = "",
                hold: bool = False,
                context: dict[str, str] = context,
                locators: dict[str, str] = locators,
            ) -> None:
                if value.casefold() in MISSING:
                    return
                item = {
                    **context,
                    "field": field,
                    "value": value,
                    "source_locator": locators[original],
                    "notes": "; ".join(filter(None, [context["notes"], note])),
                }
                (held if hold else candidates).append(item)

            mapping = (
                {
                    "collection_date": "collection_date",
                    "locality": "geo_loc_name",
                    "collection_id": "sample_title",
                    "tissue": "tissue",
                    "specimen_voucher": "specimen_voucher",
                    "life_stage": "dev_stage",
                }
                if study == "Ahuja2024"
                else {
                    "collection_date": "date_collected",
                    "locality": "location",
                    "ocean_region": "ocean",
                    "collection_id": "collection_ID",
                    "specimen_voucher": "YPM ID",
                }
            )
            for field, original in mapping.items():
                if values.get(original):
                    add(field, values[original], original)
            original = "lat_lon" if study == "Ahuja2024" else "lat_long"
            if values.get(original, "").casefold() not in MISSING:
                lat, lon = coordinates(values[original])
                add("lat_long_raw", values[original], original)
                add("latitude", lat, original)
                add("longitude", lon, original)
            if values.get("depth", "").casefold() not in MISSING:
                add("depth_original", values["depth"], "depth")
                parsed, note = depths(values["depth"])
                # These unusually deep Nanomia records need an external check;
                # preserve the source values without promoting them as metres.
                hold = study == "Ahuja2026" and values.get("sample_ID") in {
                    "WS5",
                    "WS6",
                }
                if hold:
                    note += (
                        "; pending_external_verification_of_collection_depth_and_units"
                    )
                for field, value in parsed:
                    add(field, value, "depth", note, hold=hold)
                add(
                    "collection_depth_source",
                    f"{source.relative_to(repo)}#{sheet.title}!{cells['depth']}",
                    "depth",
                    note,
                )

    candidates.sort(
        key=lambda row: (row["library_id"], row["field"], row["source_locator"])
    )
    write_tsv(out / "sheet_metadata.tsv", candidates, FIELDS)
    write_tsv(out / "sheet_metadata_pending.tsv", held, FIELDS)
    write_tsv(out / "sheet_row_matches.tsv", joins, list(joins[0]))
    (out / "sheet_raw_rows.json").write_text(
        json.dumps(raw, indent=2, ensure_ascii=False) + "\n"
    )
    (out / "sheet_sources.json").write_text(
        json.dumps(
            {
                "openpyxl_version": openpyxl.__version__,
                "sources": sources,
                "canonical_manifest_sha256": hashlib.sha256(
                    (repo / "manifest.csv").read_bytes()
                ).hexdigest(),
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n"
    )
    summary = {
        "source_rows": len(raw),
        "matched_DNA_rows": sum(row["status"] == "matched" for row in joins),
        "unique_libraries": len({row["library_id"] for row in candidates}),
        "candidates": len(candidates),
        "pending_candidates": len(held),
        "excluded_non_DNA_rows": sum(
            row["status"] == "excluded_non_DNA" for row in joins
        ),
        "unmatched_rows": [
            row
            for row in joins
            if row["status"] in {"unmatched", "conflicting_identifiers"}
        ],
    }
    (out / "sheet_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
