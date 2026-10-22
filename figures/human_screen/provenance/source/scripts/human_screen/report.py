"""Combine explicit sample summaries and evaluate the known-origin controls."""

import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--controls", action="store_true")
    parser.add_argument("summaries", nargs="+")
    args = parser.parse_args()
    rows = [json.loads(Path(p).read_text()) for p in args.summaries]
    if len({r["sample_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate sample summaries")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    if args.controls:
        human = [r for r in rows if r["sample_id"] == "control_human"]
        negative = [r for r in rows if r["sample_id"] != "control_human"]
        if len(human) != 1 or len(negative) != 3:
            raise ValueError("Expected one human and three siphonophore controls")
        sensitivity = (
            human[0]["stringent_human_compatible_pairs"] / human[0]["qc_pairs"]
        )
        failures = sum(r["stringent_human_compatible_pairs"] for r in negative)
        if sensitivity < 0.70 or failures:
            raise ValueError(
                f"Controls need review: human recovery={sensitivity:.3%}, nonhuman pairs={failures}"
            )
        args.output.with_suffix(".passed.json").write_text(
            json.dumps(
                {
                    "human_recovery": sensitivity,
                    "nonhuman_false_positives": failures,
                    "negative_pairs": sum(r["qc_pairs"] for r in negative),
                },
                indent=2,
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
