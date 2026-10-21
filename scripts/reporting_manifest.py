#!/usr/bin/env python3
"""Create a reporting-only manifest using reviewed metadata and frozen read identities."""

import argparse
import csv
import io
from pathlib import Path

from export_sample_metadata import export_metadata, read_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = export_metadata(args.output_dir)
    original = read_table(Path("config/post_expansion_manifest.csv"))
    fields = list(original[0])
    metadata = {r["library_id"]: r for r in rows}
    if set(metadata) != {r["library_id"] for r in original}:
        raise ValueError("Reporting and mapping cohorts differ")
    protected = {
        "library_id",
        "specimen_id",
        "r1_paths",
        "r2_paths",
        "read_pairs",
        "host_reference",
        "sequencing_batches",
        "include_primary",
    }
    for row in original:
        reviewed = metadata[row["library_id"]]
        for field in fields:
            if field in reviewed and field not in protected:
                row[field] = reviewed[field]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(original)
    path = args.output_dir / "manifest.csv"
    content = buffer.getvalue().encode()
    if not path.exists() or path.read_bytes() != content:
        path.write_bytes(content)


if __name__ == "__main__":
    main()
