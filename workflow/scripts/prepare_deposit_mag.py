#!/usr/bin/env python3
"""Apply reviewed sequence edits to a MAG; preserve the catalog FASTA unchanged."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA identifier: {path}: {current}")
                records[current] = []
            elif line.strip():
                if current is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                records[current].append(line.strip().upper())
    sequences = {name: "".join(parts) for name, parts in records.items()}
    for name, sequence in sequences.items():
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Empty or invalid nucleotide sequence: {path}: {name}")
    if not sequences:
        raise ValueError(f"Empty FASTA: {path}")
    return sequences


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def curate_mag(
    sequences: dict[str, str], removals: list[dict[str, str]],
    trims: list[dict[str, str]], repairs: list[dict[str, str]],
) -> tuple[dict[str, str], list[dict[str, object]]]:
    """Use original coordinates and reject ambiguous, absent or stale edit records."""
    short_ids = {name.split("__")[-1]: name for name in sequences}
    if len(short_ids) != len(sequences):
        raise ValueError("Original contig IDs are not unique within the MAG")
    edited = sequences.copy()
    events = []
    affected = set()
    for action, records in (("remove", removals), ("trim", trims), ("repair", repairs)):
        for row in records:
            contig = row["contig_id"]
            if contig not in short_ids or contig in affected:
                raise ValueError(f"Absent or multiply edited contig: {contig}")
            affected.add(contig)
            name = short_ids[contig]
            sequence = edited[name]
            original_length = len(sequence)
            if action == "remove":
                if row.get("length") and int(row["length"]) != original_length:
                    raise ValueError(f"Reviewed length changed: {contig}")
                del edited[name]
                retained_length = 0
            elif action == "trim":
                start, end = int(row["start_1based"]), int(row["end_1based"])
                if original_length != int(row["original_length"]):
                    raise ValueError(f"Reviewed length changed: {contig}")
                if not 1 <= start <= end <= original_length or not (start == 1 or end == original_length):
                    raise ValueError(f"Only reviewed terminal adaptor trimming is supported: {contig}")
                edited[name] = sequence[:start - 1] + sequence[end:]
                retained_length = len(edited[name])
            else:
                if row["action"] != "retain_first_half_of_exact_reverse_complement_palindrome":
                    raise ValueError(f"Unknown repair: {row['action']}")
                if original_length != int(row["original_length"]) or original_length % 2:
                    raise ValueError(f"Reviewed palindrome length changed: {contig}")
                if sequence != reverse_complement(sequence):
                    raise ValueError(f"Reviewed contig is not an exact palindrome: {contig}")
                edited[name] = sequence[:int(row["retained_length"])]
                retained_length = len(edited[name])
                if retained_length * 2 != original_length:
                    raise ValueError(f"Repair must retain exactly one half: {contig}")
            events.append({"contig_id": name, "action": action, "original_length": original_length,
                           "retained_length": retained_length, "reason": row["reason"]})
    if not edited or any(len(seq) < 200 for seq in edited.values()):
        raise ValueError("Curated MAG is empty or has a sequence shorter than 200 bp")
    return edited, events


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "removals", "trims", "repairs", "output", "report"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--mag", required=True)
    args = parser.parse_args()
    source = read_fasta(args.input)
    decisions = [[r for r in read_table(path) if r["mag_id"] == args.mag]
                 for path in (args.removals, args.trims, args.repairs)]
    sequences, edits = curate_mag(source, *decisions)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(".fa.part")
    with temporary.open("w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    temporary.replace(args.output)
    total = sum(map(len, sequences.values()))
    report = {"mag_id": args.mag, "input": str(args.input), "output": str(args.output),
              "input_sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
              "output_sha256": hashlib.sha256(args.output.read_bytes()).hexdigest(),
              "original_size": sum(map(len, source.values())), "genome_size": total,
              "contigs": len(sequences), "min_contig_length": min(map(len, sequences.values())),
              "gc_fraction": sum(s.count("G") + s.count("C") for s in sequences.values()) / total,
              "remaining_exact_palindromes": [name for name, s in sequences.items() if s == reverse_complement(s)],
              "edits": edits,
              "decision_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in (args.removals, args.trims, args.repairs)}}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
