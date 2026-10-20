#!/usr/bin/env python3
"""Verify local compressed FASTQ bytes against public NCBI original uploads."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.files.open() as handle:
        rows = [
            r
            for r in csv.DictReader(handle, delimiter="\t")
            if r["analysis_sample_id"] == args.sample
        ]
    assert {r["mate"] for r in rows} == {"1", "2"} and len(rows) == 2
    for row in rows:
        path = Path(row["path"])
        assert path.stat().st_size == int(row["bytes"]), path
        md5, sha = hashlib.md5(), hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
                md5.update(block)
                sha.update(block)
        if md5.hexdigest() != row["ncbi_original_md5"].zfill(32):
            raise ValueError(f"NCBI original-upload checksum mismatch: {path}")
        row["local_md5"], row["local_sha256"] = md5.hexdigest(), sha.hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps({"status": "PASS", "sample": args.sample, "files": rows}, indent=2)
        + "\n"
    )


if __name__ == "__main__":
    main()
