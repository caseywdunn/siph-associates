#!/usr/bin/env python3
"""Create depth candidates under the curator's explicit Physalia convention.

This creates a review input only. Apply it through merge_sample_metadata.py;
the accepted overlay then preserves the scoring on ordinary manifest rebuilds.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from merge_sample_metadata import CANDIDATE_COLUMNS, is_missing, read_table, write_table

CONVENTION = "docs/metadata_depth_conventions.md"
DEPTH_RECORD_FIELDS = ("depth_m", "depth_min_m", "depth_max_m", "depth_original")


def depth_candidates(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Score missing Physalia depths and distinguish records from assignments."""
    candidates = []
    for row in rows:
        has_record = any(
            not is_missing(row.get(field)) for field in DEPTH_RECORD_FIELDS
        )
        values = {}
        if row["species_current"].split()[0] == "Physalia" and not has_record:
            values = {
                "depth_m": "0",
                "collection_depth_source": CONVENTION,
                "collection_depth_basis": "curator_assigned_surface",
            }
        elif is_missing(row.get("collection_depth_basis")) and any(
            not is_missing(row.get(field)) for field in DEPTH_RECORD_FIELDS[:3]
        ):
            values = {"collection_depth_basis": "source_record"}
        for field, value in values.items():
            candidates.append(
                {
                    "library_id": row["library_id"],
                    "field": field,
                    "value": value,
                    "source_type": "curator_scoring",
                    "source": CONVENTION,
                    "source_locator": row["library_id"],
                    "source_record": "User instruction, 2026-10-04",
                    "notes": (
                        "Physalia without an explicit depth record scored 0 m by user "
                        "instruction; this does not establish beach collection."
                        if not has_record
                        else "Basis annotation only; original depth values and sources retained."
                    ),
                }
            )
    return candidates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    _, rows = read_table(args.manifest, ",")
    candidates = depth_candidates(rows)
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite review input: {args.output}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_table(args.output, CANDIDATE_COLUMNS, candidates, "\t")
    print(f"Wrote {len(candidates)} depth candidates to {args.output}")


if __name__ == "__main__":
    main()
