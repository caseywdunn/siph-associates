#!/usr/bin/env python3
"""Create a reporting-only manifest using reviewed metadata and frozen read identities."""

import argparse
import csv
import io
from pathlib import Path

from export_sample_metadata import export_metadata, read_table, unique_index


def build_reporting_manifest(
    original: list[dict[str, str]],
    reviewed: list[dict[str, str]],
    cohort: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Join biological metadata and current host routes without changing read identity."""
    metadata = unique_index(reviewed, "library_id")
    analysis_metadata = unique_index(reviewed, "analysis_sample_id")
    samples = unique_index(cohort, "sample_id")
    if set(metadata) != set(unique_index(original, "library_id")) or set(
        analysis_metadata
    ) != set(samples):
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
    rows = [dict(row) for row in original]
    for row in rows:
        sample_metadata = metadata[row["library_id"]]
        for field in row:
            if field in sample_metadata and field not in protected:
                row[field] = sample_metadata[field]
        # The original manifest predates the Halistemma reference. Processing
        # routes come from the selected cohort, never from biological metadata.
        row["host_reference"] = samples[sample_metadata["analysis_sample_id"]][
            "host_route"
        ]
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = build_reporting_manifest(
        read_table(Path("config/post_expansion_manifest.csv")),
        export_metadata(args.output_dir),
        read_table(Path("config/post_expansion_samples.tsv")),
    )
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    path = args.output_dir / "manifest.csv"
    content = buffer.getvalue().encode()
    if not path.exists() or path.read_bytes() != content:
        path.write_bytes(content)


if __name__ == "__main__":
    main()
