#!/usr/bin/env python3
"""Use annotated target COI proteins to validate barcode reading frames."""

import argparse
import csv
from pathlib import Path

from Bio import SeqIO


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("panel", type=Path)
    args = parser.parse_args()
    with (args.panel / "annotated_metadata.tsv").open() as handle:
        targets = {
            r["reference_id"]
            for r in csv.DictReader(handle, delimiter="\t")
            if r["target"] == "1"
        }
    seen = set()
    with (args.panel / "frame_anchors.faa").open("w") as handle:
        for record in SeqIO.parse(args.panel / "annotated.faa", "fasta"):
            sequence = str(record.seq)
            if record.id in targets and sequence not in seen:
                SeqIO.write(record, handle, "fasta")
                seen.add(sequence)
    if not seen:
        raise ValueError("No annotated target COI proteins for frame validation")


if __name__ == "__main__":
    main()
