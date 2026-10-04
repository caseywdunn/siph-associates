#!/usr/bin/env python3
"""Review collection-metadata additions without rerunning sequence analyses.

Candidates are long-form TSV records (CANDIDATE_COLUMNS below). The CLI writes
a proposed manifest, a durable accepted overlay, and an audit to a NEW output
directory. It never promotes these files or changes historical freeze records.
Conflicting sources require review; priority only selects the representation
of otherwise agreeing values. The same overlay reader is used by build_manifest.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from datetime import date
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
ADDITIONAL_FIELDS = [
    "depth_min_m",
    "depth_max_m",
    "depth_original",
    "collection_depth_source",
    "collection_depth_basis",
    "biosample",
    "collection_method",
    "tissue",
    "life_stage",
    "country",
    "specimen_voucher",
]
ALLOWED_FIELDS = {
    "collection_date",
    "latitude",
    "longitude",
    "lat_long_raw",
    "depth_m",
    "locality",
    "ocean_region",
    "collection_id",
    *ADDITIONAL_FIELDS,
}
NUMERIC_FIELDS = {"latitude", "longitude", "depth_m", "depth_min_m", "depth_max_m"}
COORDINATE_FIELDS = {"latitude", "longitude", "lat_long_raw"}
DEPTH_FIELDS = {
    "depth_m",
    "depth_min_m",
    "depth_max_m",
    "depth_original",
    "collection_depth_source",
    "collection_depth_basis",
}
SOURCE_PRIORITY = {"primary_sheet": 0, "biosample": 1, "gbif": 2, "curator_scoring": 3}
CANDIDATE_COLUMNS = [
    "library_id",
    "field",
    "value",
    "source_type",
    "source",
    "source_locator",
    "source_record",
    "notes",
]
PROVENANCE_COLUMNS = [
    *CANDIDATE_COLUMNS,
    "candidate_file",
    "candidate_row",
    "source_value",
]
AUDIT_COLUMNS = [*PROVENANCE_COLUMNS, "status", "existing_value", "detail"]
MISSING_VALUES = {
    "",
    "na",
    "n/a",
    "nan",
    "none",
    "null",
    "unknown",
    "missing",
    "-",
    "—",
    "not collected",
    "not provided",
    "not available",
    "not applicable",
}


def is_missing(value: Any) -> bool:
    """Recognize empty/source missing markers while retaining numeric zero."""
    text = "" if value is None else str(value).strip()
    return text.casefold() in MISSING_VALUES or text.upper().startswith("TODO")


def normalize_value(field: str, value: Any) -> str:
    """Validate candidate units/ranges and normalize only unambiguous values."""
    text = str(value).strip()
    if field in NUMERIC_FIELDS:
        try:
            number = Decimal(text)
        except InvalidOperation as error:
            raise ValueError(
                "expected a numeric value in decimal degrees/metres"
            ) from error
        if not number.is_finite():
            raise ValueError("numeric values must be finite")
        bound = 90 if field == "latitude" else 180
        if field in {"latitude", "longitude"} and abs(number) > bound:
            raise ValueError(f"{field} outside +/-{bound} degrees")
        if field.startswith("depth_") and number < 0:
            raise ValueError("collection depth must be nonnegative metres")
        return format(number.normalize(), "f") if number else "0"
    if field == "collection_date":
        if not re.fullmatch(r"\d{4}(?:-\d{2}(?:-\d{2})?)?", text):
            raise ValueError(
                "use an ISO date, year-month, or year; do not invent precision"
            )
        parts = [int(part) for part in text.split("-")]
        date(
            parts[0],
            parts[1] if len(parts) > 1 else 1,
            parts[2] if len(parts) > 2 else 1,
        )
    return text


def read_table(path: Path, delimiter: str = "\t") -> tuple[list[str], list[dict]]:
    """Read a table, rejecting duplicate headers or ragged records."""
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = list(reader.fieldnames or [])
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"missing or duplicate headers in {path}")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"ragged records in {path}")
    return fields, rows


def read_candidates(path: Path) -> list[dict[str, str]]:
    """Retain both source cell/attribute locators and candidate-file locations."""
    fields, rows = read_table(path)
    required = set(CANDIDATE_COLUMNS[:6])
    if not required.issubset(fields):
        raise ValueError(
            f"missing candidate columns in {path}: {sorted(required - set(fields))}"
        )
    unknown = set(fields) - set(PROVENANCE_COLUMNS)
    if unknown:
        raise ValueError(
            f"unrecognized columns would lose provenance in {path}: {sorted(unknown)}"
        )
    return [
        {
            **{key: row.get(key, "") for key in PROVENANCE_COLUMNS},
            "candidate_file": row.get("candidate_file") or str(path),
            "candidate_row": row.get("candidate_row") or str(index),
        }
        for index, row in enumerate(rows, 2)
    ]


def index_manifest(rows: list[dict]) -> dict[str, dict]:
    """Require exact, unique canonical IDs; never join by a partial specimen name."""
    ids = [str(row.get("library_id", "")) for row in rows]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("manifest library_id values must be nonempty and unique")
    return dict(zip(ids, rows, strict=True))


def depth_error(row: dict) -> str:
    """Check exact depths against bounds without inferring a range midpoint."""
    values = {
        field: Decimal(str(row[field]))
        for field in ("depth_m", "depth_min_m", "depth_max_m")
        if not is_missing(row.get(field, ""))
    }
    low, high = values.get("depth_min_m"), values.get("depth_max_m")
    exact = values.get("depth_m")
    if low is not None and high is not None and low > high:
        return "depth_min_m exceeds depth_max_m"
    if exact is not None:
        if low is not None and exact < low:
            return "depth_m is below depth_min_m"
        if high is not None and exact > high:
            return "depth_m is above depth_max_m"
    return ""


def apply_updates(rows: list[dict], updates: list[dict]) -> list[dict]:
    """Apply an accepted overlay during a rebuild; fail on changed source facts.

    Existing nonmissing values retain their original representation. Immutable
    columns (including read-count depth_source) are never eligible for updates.
    """
    result = [dict(row) for row in rows]
    indexed = index_manifest(result)
    touched: set[str] = set()
    for update in updates:
        library_id, field = update["library_id"], update["field"]
        if library_id not in indexed or field not in ALLOWED_FIELDS:
            raise ValueError(f"invalid accepted update: {library_id}, {field}")
        if is_missing(update["value"]):
            raise ValueError(f"missing accepted value: {library_id}, {field}")
        value = normalize_value(field, update["value"])
        current = indexed[library_id].get(field, "")
        if not is_missing(current) and normalize_value(field, current) != value:
            raise ValueError(
                f"accepted overlay conflicts with source: {library_id}, {field}"
            )
        if is_missing(current):
            indexed[library_id][field] = value
        touched.add(library_id)
    for library_id in touched:
        error = depth_error(indexed[library_id])
        if error:
            raise ValueError(f"invalid accepted range for {library_id}: {error}")
    return result


def source_key(candidate: dict) -> tuple:
    """Prefer original collection sheets only when source values already agree."""
    return (
        SOURCE_PRIORITY[candidate["source_type"]],
        candidate["source"],
        candidate["source_locator"],
        candidate["candidate_file"],
        candidate["candidate_row"],
    )


def merge_metadata(
    rows: list[dict], candidates: list[dict], existing_updates: list[dict] | None = None
) -> tuple[list, list, list]:
    """Return a fill-only manifest, accepted provenance rows, and review audit."""
    existing_updates = existing_updates or []
    baseline = apply_updates(rows, existing_updates)
    result = [dict(row) for row in baseline]
    indexed = index_manifest(result)
    original = index_manifest(baseline)
    audits: list[dict] = []
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for candidate in candidates:
        record = {key: candidate.get(key, "") for key in PROVENANCE_COLUMNS}
        library_id, field = record["library_id"], record["field"]
        record.update(status="", existing_value="", detail="")
        audits.append(record)
        if library_id not in indexed:
            record.update(
                status="unknown_library", detail="exact canonical ID required"
            )
            continue
        record["existing_value"] = str(indexed[library_id].get(field, ""))
        if field not in ALLOWED_FIELDS:
            record.update(
                status="immutable_field",
                detail="field is outside collection metadata allowlist",
            )
            continue
        if is_missing(record["value"]):
            record.update(
                status="missing_candidate", detail="source has no usable value"
            )
            continue
        if (
            record["source_type"] not in SOURCE_PRIORITY
            or not record["source"]
            or not record["source_locator"]
        ):
            record.update(
                status="invalid_candidate",
                detail="known source_type, source and source_locator required",
            )
            continue
        try:
            record["normalized_value"] = normalize_value(field, record["value"])
        except ValueError as error:
            record.update(status="invalid_candidate", detail=str(error))
            continue
        groups[library_id, field].append(record)

    for (library_id, field), records in groups.items():
        values = {record["normalized_value"] for record in records}
        current = indexed[library_id].get(field, "")
        if field == "collection_depth_source":
            # Source URLs describe evidence, not competing biological facts.
            # Preserve the original URL on each accepted row when combining them.
            if is_missing(current):
                source_values = dict.fromkeys(
                    record["normalized_value"]
                    for record in sorted(records, key=source_key)
                )
                combined = " | ".join(source_values)
                indexed[library_id][field] = combined
                for record in records:
                    record.update(
                        status="filled", detail="combined collection-depth provenance"
                    )
                    record["source_value"] = record["value"]
                    record["value"] = combined
            else:
                for record in records:
                    record.update(
                        status="existing_provenance_retained",
                        detail="nonmissing provenance text retained; additional source is in audit",
                    )
            continue
        if not is_missing(current):
            try:
                normalized_current = normalize_value(field, current)
            except ValueError:
                normalized_current = str(current)
            for record in records:
                agrees = record["normalized_value"] == normalized_current
                record.update(
                    status="agrees_with_existing"
                    if agrees
                    else "preserved_existing_conflict",
                    detail="nonmissing manifest value retained",
                )
        elif len(values) > 1:
            for record in records:
                record.update(
                    status="candidate_conflict", detail=json.dumps(sorted(values))
                )
        else:
            selected = min(records, key=source_key)
            indexed[library_id][field] = selected["normalized_value"]
            for record in records:
                record.update(status="filled", detail="all usable candidates agree")

    # Never create mixed coordinate/depth records by accepting half a conflict.
    for library_id in indexed:
        for fields in (COORDINATE_FIELDS, DEPTH_FIELDS):
            related = [
                r
                for r in audits
                if r["library_id"] == library_id and r["field"] in fields
            ]
            error = ""
            if any(
                r["field"] in NUMERIC_FIELDS
                and r["status"]
                in {
                    "candidate_conflict",
                    "preserved_existing_conflict",
                    "invalid_candidate",
                }
                for r in related
            ):
                error = (
                    "related numeric coordinate/depth field has an unresolved conflict"
                )
            if fields == DEPTH_FIELDS and not error:
                error = depth_error(indexed[library_id])
            if error:
                for record in related:
                    if record["status"] == "filled":
                        record.update(status="related_field_conflict", detail=error)
                        field = record["field"]
                        if field in original[library_id]:
                            indexed[library_id][field] = original[library_id][field]
                        else:
                            indexed[library_id].pop(field, None)

    accepted = [dict(row) for row in existing_updates]
    accepted.extend(
        {key: record.get(key, "") for key in PROVENANCE_COLUMNS}
        for record in audits
        if record["status"] in {"filled", "agrees_with_existing"}
    )
    unique = {
        tuple(row.get(key, "") for key in PROVENANCE_COLUMNS): row for row in accepted
    }
    accepted = sorted(
        unique.values(), key=lambda r: (r["library_id"], r["field"], source_key(r))
    )
    # Exercise the exact rebuild path before any proposed products are written.
    rebuilt = apply_updates(rows, accepted)
    if rebuilt != result:
        raise ValueError("accepted overlay does not reproduce the proposed manifest")
    return result, accepted, audits


def write_table(
    path: Path, fields: list[str], rows: list[dict], delimiter: str
) -> None:
    """Write a complete reviewed derivative inside its new output directory."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=REPO / "manifest.csv")
    parser.add_argument("--candidates", type=Path, nargs="+", required=True)
    parser.add_argument(
        "--updates",
        type=Path,
        default=REPO / "data/metadata/sample_metadata_updates.tsv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="new review directory; files are never promoted automatically",
    )
    parser.add_argument("--expected-libraries", type=int, default=205)
    args = parser.parse_args()
    fields, rows = read_table(args.manifest, delimiter=",")
    index_manifest(rows)
    if len(rows) != args.expected_libraries:
        parser.error(f"expected {args.expected_libraries} libraries, found {len(rows)}")
    existing = read_candidates(args.updates) if args.updates.exists() else []
    candidates = [row for path in args.candidates for row in read_candidates(path)]
    result, accepted, audit = merge_metadata(rows, candidates, existing)
    # Match the durable build_manifest schema, including fields still missing
    # throughout this cohort, so the review product and a rebuild agree.
    out_fields = fields + [field for field in ADDITIONAL_FIELDS if field not in fields]
    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_table(args.output_dir / "manifest.csv", out_fields, result, ",")
    write_table(
        args.output_dir / "sample_metadata_updates.tsv",
        PROVENANCE_COLUMNS,
        accepted,
        "\t",
    )
    write_table(
        args.output_dir / "metadata_merge_audit.tsv", AUDIT_COLUMNS, audit, "\t"
    )
    print(
        json.dumps(
            {
                "libraries": len(result),
                "candidate_statuses": Counter(r["status"] for r in audit),
                "accepted_source_records": len(accepted),
                "output_dir": str(args.output_dir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
