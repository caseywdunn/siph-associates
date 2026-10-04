#!/usr/bin/env python3
"""Apply explicit, exactly keyed review decisions to metadata candidates.

This step makes no biological decisions: each exclusion or replacement must
appear in the supplied review_decisions.tsv. Original candidates and canonical
metadata are never modified. Outputs are written only to a new review directory.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

from merge_sample_metadata import (
    PROVENANCE_COLUMNS,
    read_candidates,
    read_table,
    write_table,
)

KEY_FIELDS = ["library_id", "field", "source", "source_locator"]
DECISION_COLUMNS = [
    *KEY_FIELDS,
    "original_value",
    "action",
    "replacement_value",
    "reason",
]
AUDIT_COLUMNS = [
    *PROVENANCE_COLUMNS,
    "decision_file",
    "decision_row",
    "original_value",
    "action",
    "replacement_value",
    "reason",
]


def read_decisions(path: Path) -> list[dict[str, str]]:
    """Require the complete decision schema and retain its physical row numbers."""
    fields, rows = read_table(path)
    if set(fields) != set(DECISION_COLUMNS):
        raise ValueError(f"decision columns must be exactly {DECISION_COLUMNS}: {path}")
    return [
        {**row, "decision_file": str(path), "decision_row": str(index)}
        for index, row in enumerate(rows, 2)
    ]


def apply_decisions(
    candidates: list[dict[str, str]], decisions: list[dict[str, str]]
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Apply only unique exact matches, preserving candidate order and evidence.

    All decisions are checked against the original input. Multiple decisions
    targeting one row are rejected, rather than interpreted as sequential edits.
    """
    index: dict[tuple, list[int]] = defaultdict(list)
    for number, candidate in enumerate(candidates):
        key = tuple(candidate.get(field, "") for field in KEY_FIELDS)
        index[key].append(number)
    selected: dict[int, dict[str, str]] = {}
    for decision in decisions:
        if any(not str(decision.get(field, "")).strip() for field in KEY_FIELDS):
            raise ValueError("review decisions require complete exact source keys")
        action = decision.get("action", "")
        if action not in {"exclude", "replace"}:
            raise ValueError(f"unknown review action: {action}")
        if not str(decision.get("reason", "")).strip():
            raise ValueError("each review decision requires a reason")
        replacement = decision.get("replacement_value", "")
        if action == "replace" and not replacement.strip():
            raise ValueError("replacement_value is required for replace")
        if action == "exclude" and replacement:
            raise ValueError("exclude must have an empty replacement_value")
        key = tuple(decision[field] for field in KEY_FIELDS)
        possible = index.get(key, [])
        if not possible:
            raise ValueError(f"unmatched review decision: {key}")
        matches = [
            number
            for number in possible
            if candidates[number].get("value", "") == decision.get("original_value", "")
        ]
        if not matches:
            raise ValueError(f"stale original_value in review decision: {key}")
        if len(matches) != 1:
            raise ValueError(f"ambiguous review decision: {key}")
        number = matches[0]
        if number in selected:
            raise ValueError(f"duplicate review decisions for one candidate: {key}")
        selected[number] = decision

    reviewed, audit = [], []
    for number, candidate in enumerate(candidates):
        row = dict(candidate)
        decision = selected.get(number)
        if decision is None:
            reviewed.append(row)
            continue
        row["source_value"] = row.get("source_value") or row["value"]
        annotation = f"review {decision['action']}: {decision['reason']}"
        row["notes"] = "; ".join(
            value for value in (row.get("notes", ""), annotation) if value
        )
        if decision["action"] == "replace":
            row["value"] = decision["replacement_value"]
            reviewed.append(row)
        audit.append(
            {
                **{field: row.get(field, "") for field in PROVENANCE_COLUMNS},
                **{
                    field: decision.get(field, "")
                    for field in AUDIT_COLUMNS
                    if field not in PROVENANCE_COLUMNS
                },
            }
        )
    return reviewed, audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, nargs="+", required=True)
    parser.add_argument("--decisions", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="new directory; never promoted automatically",
    )
    args = parser.parse_args()
    candidates = [row for path in args.candidates for row in read_candidates(path)]
    decisions = read_decisions(args.decisions)
    reviewed, audit = apply_decisions(candidates, decisions)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_table(
        args.output_dir / "reviewed_candidates.tsv", PROVENANCE_COLUMNS, reviewed, "\t"
    )
    write_table(
        args.output_dir / "review_decision_audit.tsv", AUDIT_COLUMNS, audit, "\t"
    )
    print(
        json.dumps(
            {
                "input_candidates": len(candidates),
                "reviewed_candidates": len(reviewed),
                "decisions": dict(Counter(row["action"] for row in audit)),
                "output_dir": str(args.output_dir),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
