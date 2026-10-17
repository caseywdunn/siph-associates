#!/usr/bin/env python3
"""Retain collecting-locality records linked from verified LUX voucher records."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

from check_public_ypm_vouchers import catalog_number, lux_catalogs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    audit, out = args.audit_dir, args.output_dir
    (out / "raw").mkdir(parents=True, exist_ok=True)
    requests = {r["url"]: r for r in json.loads((audit / "requests.json").read_text())}
    links = []
    with (audit / "voucher_public_access.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for view in filter(None, row["lux_urls"].split(";")):
                url = view.replace("/view/", "/data/")
                body = (audit / requests[url]["file"]).read_bytes()
                if hashlib.sha256(body).hexdigest() != requests[url]["sha256"]:
                    raise ValueError(f"Source checksum differs: {url}")
                record = json.loads(body)
                if catalog_number(row["voucher_id"]) not in lux_catalogs(record):
                    raise ValueError(f"Voucher mismatch: {url}")
                for event in record.get("encountered_by", []):
                    for place in event.get("took_place_at", []):
                        links.append(
                            {
                                "library_id": row["library_id"],
                                "voucher_id": row["voucher_id"],
                                "object_url": url,
                                "locality_url": place["id"],
                            }
                        )
    receipts = []

    def fetch(url: str) -> None:
        for attempt in range(4):
            try:
                request = urllib.request.Request(
                    url,
                    headers={
                        "Accept": "application/json",
                        "User-Agent": "Mozilla/5.0 (public specimen metadata review)",
                    },
                )
                with urllib.request.urlopen(request, timeout=40) as response:
                    body = response.read()
                    resolved = response.url
                record = json.loads(body)
                if record.get("type") != "Place" or record.get("id") != url:
                    raise ValueError(f"Unexpected locality response: {url}")
                name = "raw/" + hashlib.sha256(url.encode()).hexdigest() + ".json"
                (out / name).write_bytes(body)
                receipts.append(
                    {
                        "url": url,
                        "resolved_url": resolved,
                        "file": name,
                        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                        "sha256": hashlib.sha256(body).hexdigest(),
                    }
                )
                return
            except (urllib.error.URLError, TimeoutError, ValueError) as exc:
                if attempt == 3:
                    receipts.append({"url": url, "error": str(exc)})
                    raise
                time.sleep(2**attempt)

    try:
        urls = sorted({row["locality_url"] for row in links})
        with ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(fetch, urls))
        with (out / "specimen_localities.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["library_id", "voucher_id", "object_url", "locality_url"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(links)
        print(f"Retained {len(urls)} localities for {len(links)} voucher records")
    finally:
        (out / "requests.json").write_text(
            json.dumps(sorted(receipts, key=lambda r: r["url"]), indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
