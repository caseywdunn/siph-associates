#!/usr/bin/env python3
"""Audit public GBIF and Yale LUX records for the publication's YPM IZ vouchers.

Reads specimen metadata only. Searches are candidates, not identity evidence:
every matched record must contain the exact institution/collection/catalog ID.
Raw responses, URLs, timestamps, and checksums are retained for the audit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LUX = "https://lux.collections.yale.edu"
IZ = "Invertebrate Zoology Collection, Yale Peabody Museum"


def catalog_number(value: str) -> int | None:
    """Permit known YPM IZ separators and zero padding; never match another collection."""
    match = re.fullmatch(r"YPM[ :_-]+IZ[ :_-]+0*(\d+)", value.strip(), re.I)
    return int(match[1]) if match else None


def lux_catalogs(record: dict) -> set[int]:
    return {
        number
        for identifier in record.get("identified_by", [])
        if identifier.get("type") == "Identifier"
        for number in [catalog_number(identifier.get("content", ""))]
        if number is not None
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Reuse checksummed responses from this audit.",
    )
    args = parser.parse_args()
    out = args.output_dir
    (out / "raw").mkdir(parents=True, exist_ok=True)
    source = ROOT / "figures/cohort/sample_metadata.tsv"
    with source.open() as handle:
        cohort = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if catalog_number(row["voucher_id"]) is not None
        ]
    numbers = sorted({catalog_number(row["voucher_id"]) for row in cohort})
    receipts = []
    previous = {}
    if args.resume and (out / "requests.json").exists():
        previous = {
            receipt["url"]: receipt
            for receipt in json.loads((out / "requests.json").read_text())
            if "error" not in receipt
        }

    def get(url: str) -> dict:
        filename = "raw/" + hashlib.sha256(url.encode()).hexdigest() + ".json"
        if url in previous and (out / filename).exists():
            body = (out / filename).read_bytes()
            if hashlib.sha256(body).hexdigest() == previous[url]["sha256"]:
                receipts.append(previous[url])
                return json.loads(body)
        for attempt in range(4):
            try:
                request = urllib.request.Request(
                    url,
                    headers={
                        "User-Agent": "Mozilla/5.0 (compatible; luxy/0.0.7; public specimen metadata audit)",
                        "Accept": "application/json",
                    },
                )
                with urllib.request.urlopen(request, timeout=40) as response:
                    body = response.read()
                    resolved = response.url
                data = json.loads(body)
                if not isinstance(data, dict):
                    raise ValueError(f"Expected JSON object: {url}")
                (out / filename).write_bytes(body)
                receipts.append(
                    {
                        "url": url,
                        "resolved_url": resolved,
                        "file": filename,
                        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                        "sha256": hashlib.sha256(body).hexdigest(),
                        "attempt": attempt + 1,
                    }
                )
                return data
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                if attempt == 3:
                    receipts.append({"url": url, "error": str(exc)})
                    raise
                time.sleep(2**attempt)
        raise RuntimeError("Unreachable retrieval state")

    def lux_search(
        batch: list[int], collection_scoped: bool = True
    ) -> tuple[list[dict], list[str]]:
        values = sorted({value for n in batch for value in (str(n), f"{n:06d}")})
        text_query = {"OR": [{"text": v} for v in values]}
        query = (
            {"AND": [{"memberOf": {"name": IZ}}, text_query]}
            if collection_scoped
            else text_query
        )
        url = (
            LUX + "/api/search/item?" + urllib.parse.urlencode({"q": json.dumps(query)})
        )
        items, urls = [], []
        while url:
            if url in urls:
                raise ValueError("Repeated LUX pagination URL")
            urls.append(url)
            data = get(url)
            items.extend(data["orderedItems"])
            url = data.get("next", {}).get("id", "")
        return items, urls

    def gbif_search(batch: list[int]) -> tuple[list[dict], list[str]]:
        params = [
            ("institutionCode", "YPM"),
            ("collectionCode", "IZ"),
            ("limit", "300"),
        ]
        for number in batch:
            variants = {str(number), f"{number:06d}"}
            for value in sorted(variants):
                for catalog in [
                    f"YPM IZ {value}",
                    f"YPM:IZ:{value}",
                    f"YPM-IZ-{value}",
                    f"IZ.{value}",
                    value,
                ]:
                    params.append(("catalogNumber", catalog))
        offset, items, urls = 0, [], []
        while True:
            url = "https://api.gbif.org/v1/occurrence/search?" + urllib.parse.urlencode(
                [*params, ("offset", str(offset))]
            )
            urls.append(url)
            data = get(url)
            items.extend(data["results"])
            if data["endOfRecords"]:
                return items, urls
            if not data["results"]:
                raise ValueError("GBIF pagination stopped before endOfRecords")
            offset += len(data["results"])

    def gbif_catalog(record: dict) -> int | None:
        number = catalog_number(record.get("catalogNumber", ""))
        if (
            number is None
            and record.get("institutionCode") == "YPM"
            and record.get("collectionCode") == "IZ"
        ):
            match = re.fullmatch(r"(?:IZ\.)?0*(\d+)", record.get("catalogNumber", ""))
            number = int(match[1]) if match else None
        return number

    try:
        config = get(LUX + "/api/advanced-search-config")
        assert {"memberOf", "text"} <= config["terms"]["item"].keys()
        # Controls validate the complete search -> record -> catalog-ID route.
        control_items, _ = lux_search([35347])
        control_records = [get(item["id"]) for item in control_items]
        assert any(35347 in lux_catalogs(record) for record in control_records)
        gbif_control, _ = gbif_search([116025])
        assert any(gbif_catalog(record) == 116025 for record in gbif_control)
        print("Both public-record controls passed", flush=True)

        matches = {"gbif": {n: [] for n in numbers}, "lux": {n: [] for n in numbers}}
        queries = {"gbif": {}, "lux": {}}
        jobs = [
            (provider, numbers[start : start + 14])
            for provider in ("gbif", "lux")
            for start in range(0, len(numbers), 14)
        ]

        def search(job):
            provider, batch = job
            items, urls = (gbif_search if provider == "gbif" else lux_search)(batch)
            print(
                f"{provider}: searched {len(batch)} vouchers; {len(items)} candidates",
                flush=True,
            )
            return provider, batch, items, urls

        lux_items = {}
        with ThreadPoolExecutor(max_workers=3) as pool:
            for provider, batch, items, urls in pool.map(search, jobs):
                for number in batch:
                    queries[provider][number] = list(urls)
                if provider == "lux":
                    for item in items:
                        lux_items[item["id"]] = item
                else:
                    for item in items:
                        number = gbif_catalog(item)
                        if number in matches[provider]:
                            matches[provider][number].append(item)
        print(f"Retrieving {len(lux_items)} distinct LUX candidate records", flush=True)
        with ThreadPoolExecutor(max_workers=3) as pool:
            for index, record in enumerate(pool.map(get, sorted(lux_items)), 1):
                for number in lux_catalogs(record) & set(numbers):
                    matches["lux"][number].append(record)
                if index % 20 == 0:
                    print(
                        f"Verified catalog identifiers in {index} LUX records",
                        flush=True,
                    )

        # Search all LUX collections before concluding that a voucher has no match.
        # This catches records with absent or unexpected collection membership.
        missing = [n for n in numbers if not matches["lux"][n]]
        fallback_items = {}
        for start in range(0, len(missing), 14):
            batch = missing[start : start + 14]
            items, urls = lux_search(batch, collection_scoped=False)
            for number in batch:
                queries["lux"][number].extend(urls)
            for item in items:
                if item["id"] not in lux_items:
                    fallback_items[item["id"]] = item
        print(
            f"Broader LUX search: {len(missing)} unmatched vouchers; "
            f"{len(fallback_items)} additional candidate records",
            flush=True,
        )
        with ThreadPoolExecutor(max_workers=3) as pool:
            for record in pool.map(get, sorted(fallback_items)):
                for number in lux_catalogs(record) & set(numbers):
                    matches["lux"][number].append(record)

        fields = [
            "library_id",
            "voucher_id",
            "species_current",
            "gbif_status",
            "gbif_urls",
            "lux_status",
            "lux_urls",
            "lux_record_labels",
            "gbif_search_urls",
            "lux_search_urls",
        ]
        rows = []
        for specimen in cohort:
            number = catalog_number(specimen["voucher_id"])
            gbif, lux = matches["gbif"][number], matches["lux"][number]
            rows.append(
                {
                    **{key: specimen[key] for key in fields[:3]},
                    "gbif_status": "public_record_found" if gbif else "no_exact_match",
                    "gbif_urls": ";".join(
                        sorted(
                            {
                                f"https://www.gbif.org/occurrence/{r['key']}"
                                for r in gbif
                            }
                        )
                    ),
                    "lux_status": "public_record_found" if lux else "no_exact_match",
                    "lux_urls": ";".join(
                        sorted({r["id"].replace("/data/", "/view/") for r in lux})
                    ),
                    "lux_record_labels": " | ".join(r.get("_label", "") for r in lux),
                    "gbif_search_urls": " | ".join(queries["gbif"][number]),
                    "lux_search_urls": " | ".join(queries["lux"][number]),
                }
            )
        with (out / "voucher_public_access.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        with (out / "unconfirmed_vouchers.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fields[:3],
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()
            writer.writerows(
                row
                for row in rows
                if row["gbif_status"] == row["lux_status"] == "no_exact_match"
            )
        summary = {
            "checked_utc": datetime.now(timezone.utc).isoformat(),
            "libraries": len(cohort),
            "unique_ypm_vouchers": len(numbers),
            "gbif_public": sum(bool(v) for v in matches["gbif"].values()),
            "lux_public": sum(bool(v) for v in matches["lux"].values()),
            "not_found_in_either": [
                f"YPM:IZ:{n}"
                for n in numbers
                if not matches["lux"][n] and not matches["gbif"][n]
            ],
            "controls": {"lux_YPM_IZ_035347": "passed", "gbif_YPM_IZ_116025": "passed"},
            "input_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "method": "LUX collection-scoped full-text catalog-number discovery, both zero-padded and unpadded, all pages, followed by exact Identifier validation; unmatched vouchers also searched across all LUX collections. GBIF institution/collection-scoped catalog-number variants, all pages, exact catalog validation.",
            "limitation": "No exact match is a search result, not a claim that the museum record is private or does not exist.",
        }
        (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2), flush=True)
    finally:
        (out / "requests.json").write_text(
            json.dumps(sorted(receipts, key=lambda r: r["url"]), indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
