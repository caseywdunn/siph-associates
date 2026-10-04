"""Normalize retained BioSample XML without modifying the project manifest.

Run from the analysis repository root with Python 3.10 or newer. Network
retrieval is separate: biosample_requests.json records the exact batch URLs.
Only depth values carrying explicit metre units become numeric candidates.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET


BASE = Path(__file__).resolve().parent
MISSING = {"", "missing", "not applicable", "not collected", "not provided", "na", "n/a"}
COLS = [
    "library_id", "field", "value", "source_type", "source", "source_locator",
    "source_record", "notes",
]
METRES = r"(?:m|meters?|metres?)"
NUMBER = r"\d+(?:\.\d+)?"


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    """Write a deterministic tabular derivative of the source XML."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def numeric(value: str) -> str:
    """Represent a source decimal without changing its numerical value."""
    return format(float(value), ".12g")


def main() -> None:
    """Validate retrieved records and emit candidates with source provenance."""
    mapping = list(csv.DictReader((BASE / "biosample_mapping.tsv").open(), delimiter="\t"))
    requests = json.loads((BASE / "biosample_requests.json").read_text())
    records: dict[str, ET.Element] = {}
    source_paths: dict[str, str] = {}
    receipts = []
    for request in requests:
        path = BASE / f"{request['batch']}.xml"
        data = path.read_bytes()
        samples = ET.fromstring(data).findall(".//BioSample")
        found = {sample.attrib["accession"] for sample in samples}
        if found != set(request["accessions"]):
            raise ValueError(f"Incomplete or unexpected BioSample batch: {path}")
        for sample in samples:
            accession = sample.attrib["accession"]
            if accession in records:
                raise ValueError(f"Duplicate BioSample accession: {accession}")
            records[accession] = sample
            source_paths[accession] = path.name
        receipts.append({
            "path": path.name,
            "url": request["url"],
            "retrieved_utc_from_mtime": datetime.fromtimestamp(
                path.stat().st_mtime, timezone.utc
            ).isoformat(),
            "sha256": hashlib.sha256(data).hexdigest(),
            "bytes": len(data),
            "biosample_records": len(samples),
        })

    candidates: list[dict[str, str]] = []
    attributes: list[dict[str, str]] = []
    issues: list[dict[str, str]] = []
    coverage: list[dict[str, str]] = []
    for library in mapping:
        library_id = library["library_id"]
        accession = library["biosample"]
        if not accession:
            coverage.append({"library_id": library_id, "biosample": "", "status": library["mapping_status"]})
            continue
        sample = records[accession]
        url = f"https://www.ncbi.nlm.nih.gov/biosample/{accession}"
        coverage.append({"library_id": library_id, "biosample": accession, "status": "retrieved_exact_accession"})

        def add(field: str, value: str, locator: str, notes: str = "") -> None:
            candidates.append({
                "library_id": library_id, "field": field, "value": value,
                "source_type": "biosample", "source": url,
                "source_locator": locator, "source_record": accession,
                "notes": notes,
            })

        add("biosample", accession, "BioSample/@accession")
        for attribute in sample.findall("./Attributes/Attribute"):
            name = attribute.get("harmonized_name", attribute.get("attribute_name", ""))
            value = (attribute.text or "").strip()
            locator = f"Attributes/Attribute[@attribute_name='{attribute.get('attribute_name', '')}']"
            attributes.append({
                "library_id": library_id, "biosample": accession,
                "attribute_name": attribute.get("attribute_name", ""),
                "harmonized_name": name, "value": value, "source": url,
                "raw_response": source_paths[accession],
            })
            if value.lower() in MISSING:
                continue
            if name in {"collection_date", "tissue", "specimen_voucher"}:
                add(name, value, locator)
            elif name == "dev_stage":
                add("life_stage", value, locator)
            elif name == "geo_loc_name":
                place = value.partition(":")[0].strip()
                add("locality", value, locator, "Original geographic-location label retained; not an inference-region assignment.")
                if not place.endswith(("Ocean", "Sea")):
                    add("country", place, locator, "Top-level geographic jurisdiction retained as submitted, including named territories.")
            elif name == "lat_lon":
                add("lat_long_raw", value, locator)
                match = re.fullmatch(rf"({NUMBER})\s*([NS])\s+({NUMBER})\s*([EW])", value)
                if not match:
                    issues.append({"library_id": library_id, "field": "lat_lon", "value": value, "issue": "unparsed_coordinate_format", "source": url})
                    continue
                lat, ns, lon, ew = match.groups()
                latitude = float(lat) * (-1 if ns == "S" else 1)
                longitude = float(lon) * (-1 if ew == "W" else 1)
                if abs(latitude) > 90 or abs(longitude) > 180:
                    raise ValueError(f"Invalid geographic coordinates for {library_id}")
                add("latitude", numeric(str(latitude)), locator, "Hemisphere sign applied; no geographic correction inferred.")
                add("longitude", numeric(str(longitude)), locator, "Hemisphere sign applied; no geographic correction inferred.")
                if accession == "SAMN32803051":
                    issues.append({"library_id": library_id, "field": "longitude", "value": value, "issue": "west_longitude_conflicts_with_Guam_label_and_existing_east_coordinate", "source": url})
            elif name == "depth":
                add("depth_original", value, locator)
                add("collection_depth_source", url, locator)
                scalar = re.fullmatch(rf"({NUMBER})\s*{METRES}", value, re.I)
                interval = re.fullmatch(rf"({NUMBER})\s*[-–]\s*({NUMBER})\s*{METRES}", value, re.I)
                if scalar:
                    add("depth_m", numeric(scalar.group(1)), locator, "Explicit metre unit retained; source depth definition otherwise unchanged.")
                elif interval:
                    low, high = interval.groups()
                    if float(low) > float(high):
                        raise ValueError(f"Reversed depth interval for {library_id}")
                    add("depth_min_m", numeric(low), locator, "Explicit metre interval; no midpoint imputed.")
                    add("depth_max_m", numeric(high), locator, "Explicit metre interval; no midpoint imputed.")
                else:
                    issues.append({"library_id": library_id, "field": "depth", "value": value, "issue": "numeric_depth_with_missing_units_or_unresolved_format; no_depth_m_candidate", "source": url})
                if accession == "SAMN32803049":
                    issues.append({"library_id": library_id, "field": "depth_m", "value": "297", "issue": "BioSample_matches_2024_sheet_297m_but_2026_sheet_reports_296; retain_conflict", "source": url})

    candidates.sort(key=lambda row: (row["library_id"], row["field"], row["value"]))
    write_tsv(BASE / "biosample_candidates.tsv", candidates, COLS)
    write_tsv(BASE / "biosample_attributes.tsv", attributes, ["library_id", "biosample", "attribute_name", "harmonized_name", "value", "source", "raw_response"])
    write_tsv(BASE / "biosample_issues.tsv", issues, ["library_id", "field", "value", "issue", "source"])
    write_tsv(BASE / "biosample_coverage.tsv", coverage, ["library_id", "biosample", "status"])
    summary = {
        "normalized_utc": datetime.now(timezone.utc).isoformat(),
        "normalizer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "mapping_sha256": hashlib.sha256((BASE / "biosample_mapping.tsv").read_bytes()).hexdigest(),
        "retrievals": receipts,
        "library_count": len(mapping), "retrieved_biosamples": len(records),
        "coverage": dict(Counter(row["status"] for row in coverage)),
        "candidate_field_counts": dict(Counter(row["field"] for row in candidates)),
        "issues": len(issues),
    }
    (BASE / "biosample_provenance.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "retrievals"}, indent=2))


if __name__ == "__main__":
    main()
