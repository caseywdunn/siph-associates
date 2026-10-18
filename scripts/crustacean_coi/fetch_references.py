#!/usr/bin/env python3
"""Download public COI reference inputs and preserve their retrieval identities."""

import argparse
import hashlib
import json
import shutil
import urllib.request
from datetime import UTC, datetime
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    receipt = args.output / "downloads.json"
    previous = json.loads(receipt.read_text()) if receipt.exists() else []
    by_name = {r["filename"]: r for r in previous}
    rows = []
    for source in json.loads(args.sources.read_text()):
        dest = args.output / source["filename"]
        if dest.exists() and source["filename"] in by_name:
            old = by_name[source["filename"]]
            with dest.open("rb") as handle:
                sha = hashlib.file_digest(handle, "sha256").hexdigest()
            if old["url"] != source["url"] or old["sha256"] != sha:
                raise ValueError(
                    f"Existing download does not match its receipt: {dest}"
                )
            rows.append(old)
            continue
        # A pre-existing inspection download has no receipt and is replaced by
        # a verified complete download; partial transfers never become inputs.
        temporary = dest.with_suffix(dest.suffix + ".part")
        request = urllib.request.Request(
            source["url"], headers={"User-Agent": "siph-associates-coi-pilot"}
        )
        with (
            urllib.request.urlopen(request, timeout=120) as response,
            temporary.open("wb") as output,
        ):
            headers = dict(response.headers)
            shutil.copyfileobj(response, output)
        with temporary.open("rb") as handle:
            sha = hashlib.file_digest(handle, "sha256").hexdigest()
        if not temporary.stat().st_size:
            raise ValueError(f"Empty reference download: {source['url']}")
        temporary.replace(dest)
        rows.append(
            {
                **source,
                "retrieved_utc": datetime.now(UTC).isoformat(),
                "bytes": dest.stat().st_size,
                "sha256": sha,
                "headers": headers,
            }
        )
        receipt.write_text(json.dumps(rows, indent=2) + "\n")
        print(source["name"], dest.stat().st_size, flush=True)
    receipt.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
