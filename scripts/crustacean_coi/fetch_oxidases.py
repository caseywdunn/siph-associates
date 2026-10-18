#!/usr/bin/env python3
"""Save public prokaryotic oxidase sequences, annotations and release identity."""

import argparse
import csv
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    root = Path(config["run"]) / "source"
    root.mkdir(parents=True, exist_ok=True)
    destination = root / "uniprot_oxidases.tsv"
    receipt = root / "uniprot_receipt.json"
    if destination.exists() and receipt.exists():
        expected = json.loads(receipt.read_text())["sha256"]
        if hashlib.sha256(destination.read_bytes()).hexdigest() != expected:
            raise ValueError("Existing UniProt download does not match its receipt")
        return
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        dict(
            query=config["uniprot_query"],
            format="tsv",
            fields="accession,id,protein_name,organism_name,organism_id,lineage,reviewed,sequence",
        )
    )
    partial = destination.with_suffix(".partial")
    with (
        urllib.request.urlopen(url, timeout=300) as response,
        partial.open("wb") as out,
    ):
        headers = dict(response.headers)
        while chunk := response.read(1024 * 1024):
            out.write(chunk)
    with partial.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) < 1000 or not all(r.get("Sequence") and r.get("Entry") for r in rows):
        raise ValueError("Incomplete or unexpected UniProt oxidase response")
    if len({r["Entry"] for r in rows}) != len(rows):
        raise ValueError("Duplicate UniProt accessions")
    partial.replace(destination)
    receipt.write_text(
        json.dumps(
            dict(
                url=url,
                query=config["uniprot_query"],
                retrieved_utc=datetime.now(UTC).isoformat(),
                headers=headers,
                records=len(rows),
                sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
            ),
            indent=2,
        )
        + "\n"
    )
    print(f"Saved {len(rows):,} oxidase proteins to {destination}")


if __name__ == "__main__":
    main()
