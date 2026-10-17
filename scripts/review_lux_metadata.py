#!/usr/bin/env python3
"""Compare verified museum records with publication metadata and propose curation.

Reads retained JSON only. It never edits the curation overlay or analysis inputs.
Existing measurements and identifiers are retained when museum sources disagree.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from decimal import Decimal
from pathlib import Path

from check_public_ypm_vouchers import catalog_number, lux_catalogs
from export_sample_metadata import ROOT, missing, read_table, write_table

COUNTRIES = {
    "USA",
    "Australia",
    "France",
    "Mexico",
    "Japan",
    "Guam",
    "Bermuda",
    "New Zealand",
}


def labels(node: dict) -> set[str]:
    return {c.get("_label", "") for c in node.get("classified_as", [])}


def values(record: dict, classification: str) -> list[str]:
    return [
        n["content"]
        for n in record.get("identified_by", [])
        if classification in labels(n)
    ]


def depth_interval(text: str) -> tuple[Decimal, Decimal] | None:
    match = re.fullmatch(
        r"(\d+(?:\.\d+)?)(?:\s*[-–]\s*(\d+(?:\.\d+)?))?\s+meters", text
    )
    if not match:
        return None
    low, high = Decimal(match[1]), Decimal(match[2] or match[1])
    if low > high:
        raise ValueError(f"Reversed museum depth interval: {text}")
    return low, high


def museum_date(label: str) -> str:
    dates = re.findall(
        r"\b\d{1,2} (?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec) \d{4}\b", label
    )
    return (
        datetime.strptime(dates[0], "%d %b %Y").date().isoformat()
        if len(dates) == 1
        else ""
    )


def coordinate_status(wkt: str, row: dict) -> str:
    # A locality uncertainty polygon is not a measured point. Check its envelope;
    # never replace study coordinates with the centre of a museum geometry.
    if (
        not wkt.startswith("POLYGON((")
        or missing(row["latitude"])
        or missing(row["longitude"])
    ):
        return "not_comparable"
    pairs = re.findall(r"(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)", wkt)
    if not pairs:
        raise ValueError(f"Unparsed locality geometry: {wkt}")
    xs, ys = zip(*[(Decimal(x), Decimal(y)) for x, y in pairs])
    if max(xs) - min(xs) > 180:
        return "not_comparable"
    x, y = Decimal(row["longitude"]), Decimal(row["latitude"])
    tolerance = Decimal("0.00001")
    within = (
        min(xs) - tolerance <= x <= max(xs) + tolerance
        and min(ys) - tolerance <= y <= max(ys) + tolerance
    )
    return "within_locality_envelope" if within else "outside_locality_envelope"


def load_responses(directory: Path) -> dict[str, dict]:
    result = {}
    for receipt in json.loads((directory / "requests.json").read_text()):
        if "error" in receipt:
            raise ValueError(f"Failed source retrieval: {receipt['url']}")
        body = (directory / receipt["file"]).read_bytes()
        if hashlib.sha256(body).hexdigest() != receipt["sha256"]:
            raise ValueError(f"Source checksum differs: {receipt['file']}")
        result[receipt["url"]] = json.loads(body)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--localities-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--overlay", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    meta = {r["library_id"]: r for r in read_table(args.metadata)}
    baseline = {r["library_id"]: r for r in read_table(ROOT / "manifest.csv")}
    old_notes = defaultdict(set)
    for row in read_table(args.overlay):
        if row["field"] == "notes":
            old_notes[row["library_id"]].add(row["value"])
    objects = load_responses(args.audit_dir)
    places = load_responses(args.localities_dir)
    comparisons, proposals, snapshots = [], [], []
    notes = defaultdict(list)
    note_sources = defaultdict(set)

    def compare(
        library: str, field: str, candidate: str, status: str, source: str, locator: str
    ) -> None:
        row = meta[library]
        current = row.get(field, "")
        if field == "coordinates":
            current = f"latitude={row['latitude']}; longitude={row['longitude']}"
        elif field == "depth_m" and all(
            not missing(row[f]) for f in ["depth_min_m", "depth_max_m"]
        ):
            current = f"{current} m; reported range {row['depth_min_m']}–{row['depth_max_m']} m"
        comparisons.append(
            {
                "library_id": library,
                "field": field,
                "current_value": current,
                "museum_value": candidate,
                "status": status,
                "source": source,
                "source_locator": locator,
            }
        )

    def propose(
        library: str, field: str, value: str, source: str, locator: str, note: str
    ) -> None:
        original = baseline[library].get(field, "")
        if field == "voucher_id":
            original = baseline[library].get("specimen_voucher", "")
        elif field == "location":
            original = baseline[library].get("locality", "")
        elif field in {"notes", "sample_id", "collection_site"}:
            original = ""
        proposals.append(
            {
                "library_id": library,
                "field": field,
                "original_value": original,
                "value": value,
                "source": source,
                "source_locator": locator,
                "notes": note,
            }
        )

    def explain(library: str, text: str, source: str) -> None:
        notes[library].append(text)
        note_sources[library].add(source)

    links = read_table(args.localities_dir / "specimen_localities.tsv")
    if len({r["library_id"] for r in links}) != len(links):
        raise ValueError("Multiple museum localities require individual review")
    for link in links:
        library = link["library_id"]
        row = meta[library]
        object_url, place_url = link["object_url"], link["locality_url"]
        obj, place = objects[object_url], places[place_url]
        if catalog_number(row["voucher_id"]) not in lux_catalogs(obj):
            raise ValueError(f"Voucher mismatch: {library}")
        if place_url not in {
            p["id"]
            for e in obj.get("encountered_by", [])
            for p in e.get("took_place_at", [])
        }:
            raise ValueError(
                f"Current specimen no longer links to this locality: {library}"
            )
        compare(
            library,
            "voucher_id",
            row["voucher_id"],
            "confirmed",
            object_url,
            "identified_by;YPM External Catalog Number",
        )
        propose(
            library,
            "voucher_id",
            row["voucher_id"],
            object_url,
            "identified_by;YPM External Catalog Number",
            "Exact museum catalog identifier independently confirms the study voucher.",
        )
        name = place.get("_label", "")
        primary = values(place, "Primary Name")
        site = primary[0] if len(primary) == 1 else ""
        date = museum_date(obj.get("_label", ""))
        field_ids = [
            re.sub(r" \(field number\)$", "", v) for v in values(obj, "field number")
        ]
        personal_ids = [
            re.sub(r" \(personal specimen number\)$", "", v)
            for v in values(obj, "personal specimen number")
        ]
        other_localities = values(place, "Other Locality Number")
        depth_texts = [
            n["content"]
            for n in place.get("referred_to_by", [])
            if "Depth below sea level" in labels(n)
        ]
        geometry = place.get("defined_by", "")
        geo_status = coordinate_status(geometry, row)
        snapshots.append(
            {
                "library_id": library,
                "voucher_id": row["voucher_id"],
                "object_url": object_url,
                "locality_url": place_url,
                "object_label": obj.get("_label", ""),
                "locality_label": name,
                "primary_locality_name": site,
                "museum_field_numbers": " | ".join(field_ids),
                "museum_personal_numbers": " | ".join(personal_ids),
                "other_locality_numbers": " | ".join(other_localities),
                "depth_text": " | ".join(depth_texts),
                "geometry": geometry,
                "coordinate_status": geo_status,
                "collection_date": date,
            }
        )
        compare(library, "coordinates", geometry, geo_status, place_url, "defined_by")
        if geo_status == "outside_locality_envelope":
            explain(
                library,
                "Study coordinates fall outside the museum locality's mapped extent; study coordinates retained pending review.",
                place_url,
            )
        if date:
            status = "agrees" if date == row["collection_date"] else "disagrees"
            compare(
                library,
                "collection_date",
                date,
                status,
                object_url,
                "_label;explicit day month year",
            )
            if status == "disagrees":
                explain(
                    library,
                    f"Collection date differs: study {row['collection_date']}; museum {date}. Study date retained pending review.",
                    object_url,
                )

        # An identifier named 'field number' may identify a specimen rather than
        # an event. Use it to corroborate an existing label, never fill a blank.
        for field, candidates in [
            ("collection_id", sorted(set(field_ids + other_localities))),
            ("sample_id", personal_ids),
        ]:
            if candidates:
                status = (
                    "label_matches" if row[field] in candidates else "other_identifier"
                )
                if (
                    field == "collection_id"
                    and library.startswith("Ahuja2024:")
                    and field_ids
                    and row[field] not in field_ids
                ):
                    status = "alternate_field_label"
                compare(
                    library,
                    field,
                    " | ".join(candidates),
                    status,
                    object_url,
                    "identified_by;museum field/personal identifier;linked locality",
                )
                if status == "label_matches":
                    propose(
                        library,
                        field,
                        row[field],
                        object_url,
                        "identified_by;museum identifier;linked locality",
                        "Museum identifier cross-reference agrees with the study label; original study source establishes its sample/event meaning.",
                    )
        if (
            library.startswith("Ahuja2024:")
            and field_ids
            and row["collection_id"] not in field_ids
        ):
            explain(
                library,
                f"Museum field label is {' / '.join(field_ids)}; study collection label {row['collection_id']} retained.",
                object_url,
            )

        country_candidates = [
            c
            for c in COUNTRIES
            if re.search(r"(?:^|\. )" + re.escape(c) + r"(?:\.|$)", name)
        ]
        country = country_candidates[0] if len(country_candidates) == 1 else ""
        country_conflict = "Gulf of California" in name and country == "USA"
        if country_conflict:
            compare(
                library,
                "country",
                country,
                "internally_inconsistent_locality",
                place_url,
                "_label",
            )
            explain(
                library,
                "Museum locality lists both Gulf of California and USA; country left unresolved.",
                place_url,
            )
        elif country:
            status = (
                "agrees"
                if country == row["country"]
                else "different_label"
                if not missing(row["country"])
                else "candidate_addition"
            )
            compare(
                library,
                "country",
                country,
                status,
                place_url,
                "_label;explicit country/territory",
            )
            if missing(row["country"]) and geo_status != "outside_locality_envelope":
                propose(
                    library,
                    "country",
                    country,
                    place_url,
                    "_label;explicit country/territory",
                    "Country/territory explicitly named by the museum; not inferred from coordinates.",
                )

        # Method-only strings do not establish a geographic collection site.
        useful_site = (
            site and site != "Blue water dive" and not site.startswith("[D1398-SS8]")
        )
        compare(
            library,
            "collection_site",
            site,
            "candidate_addition"
            if missing(row["collection_site"]) and useful_site
            else "retained_source_description",
            place_url,
            "identified_by;Primary Name",
        )
        if (
            missing(row["collection_site"])
            and useful_site
            and geo_status != "outside_locality_envelope"
        ):
            propose(
                library,
                "collection_site",
                site,
                place_url,
                "identified_by;Primary Name",
                "Museum collection-locality name; original spelling retained.",
            )
            location = (
                site
                if row["location"].casefold() in site.casefold()
                else row["location"] + ": " + site
            )
            propose(
                library,
                "location",
                location,
                place_url,
                "identified_by;Primary Name",
                "Existing broad location enriched with the museum locality description.",
            )

        method = ""
        label = obj.get("_label", "")
        if "Blue water dive" in name:
            method = "blue-water diving"
        for vehicle in ["Tiburon", "Ventana", "Doc Ricketts"]:
            if f"ROV {vehicle}" in label:
                method = f"ROV {vehicle}"
        if "Johnson Sea Link II" in label:
            method = "collection using Johnson Sea Link II"
        if method and missing(row["collection_method"]):
            propose(
                library,
                "collection_method",
                method,
                object_url,
                "_label;explicit collecting platform/method",
                "Explicit museum method/platform; vehicle retained without assuming vehicle type for Johnson Sea Link II.",
            )
        for text in depth_texts:
            interval = depth_interval(text)
            if interval is None:
                compare(
                    library,
                    "depth_m",
                    text,
                    "units_or_format_unresolved",
                    place_url,
                    "referred_to_by;Depth below sea level",
                )
                continue
            low, high = interval
            point = Decimal(row["depth_m"]) if not missing(row["depth_m"]) else None
            bounds = all(not missing(row[f]) for f in ["depth_min_m", "depth_max_m"])
            conflict = (point is not None and not low <= point <= high) or (
                bounds
                and (Decimal(row["depth_min_m"]), Decimal(row["depth_max_m"]))
                != interval
            )
            compare(
                library,
                "depth_m",
                text,
                "disagrees" if conflict else "compatible",
                place_url,
                "referred_to_by;Depth below sea level",
            )
            if conflict:
                existing = (
                    f"{row['depth_min_m']}–{row['depth_max_m']} m"
                    if bounds
                    else f"{row['depth_m']} m"
                )
                explain(
                    library,
                    f"Depth differs: study {existing}; museum {text}. Study depth retained pending review.",
                    place_url,
                )
            elif (
                point is not None
                and not bounds
                and all(missing(row[f]) for f in ["depth_min_m", "depth_max_m"])
            ):
                for field, value in [("depth_min_m", low), ("depth_max_m", high)]:
                    propose(
                        library,
                        field,
                        str(value),
                        place_url,
                        "referred_to_by;Depth below sea level",
                        "Museum depth bounds added; compatible explicit study point depth retained.",
                    )
                if low != high:
                    explain(
                        library,
                        f"Museum gives a {low}–{high} m depth interval; the study's explicit {point} m point depth is retained.",
                        place_url,
                    )

    for library, descriptions in notes.items():
        value = " ".join([*sorted(old_notes[library]), *dict.fromkeys(descriptions)])
        propose(
            library,
            "notes",
            value,
            " | ".join(sorted(note_sources[library])),
            "Museum comparison; prior sourced notes retained from the input overlay",
            "Human-facing source disagreements and alternate identifiers; existing values retained.",
        )
    write_table(args.output_dir / "museum_metadata.tsv", snapshots, list(snapshots[0]))
    fields = [
        "library_id",
        "field",
        "current_value",
        "museum_value",
        "status",
        "source",
        "source_locator",
    ]
    write_table(args.output_dir / "field_comparisons.tsv", comparisons, fields)
    review_status = {
        "disagrees",
        "outside_locality_envelope",
        "internally_inconsistent_locality",
        "units_or_format_unresolved",
        "alternate_field_label",
    }
    write_table(
        args.output_dir / "review_needed.tsv",
        [r for r in comparisons if r["status"] in review_status],
        fields,
    )
    write_table(
        args.output_dir / "proposed_updates.tsv",
        proposals,
        [
            "library_id",
            "field",
            "original_value",
            "value",
            "source",
            "source_locator",
            "notes",
        ],
    )
    summary = {
        "vouchers_reviewed": len(snapshots),
        "comparison_status": dict(Counter(r["status"] for r in comparisons)),
        "proposals_by_field": dict(Counter(r["field"] for r in proposals)),
        "inputs_sha256": {
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                args.metadata,
                args.overlay,
                args.audit_dir / "requests.json",
                args.localities_dir / "requests.json",
                Path(__file__),
            ]
        },
    }
    (args.output_dir / "review_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
