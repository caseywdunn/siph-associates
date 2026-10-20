#!/usr/bin/env python3
"""Check raw-read identity, paired pilot accounting and host-removal usefulness."""

import argparse
import csv
import gzip
import json
import subprocess
from pathlib import Path


def paired_count(r1: Path, r2: Path) -> int:
    """Require complete records, matching names and sequence/quality lengths."""
    count = 0
    with gzip.open(r1, "rt") as first, gzip.open(r2, "rt") as second:
        while True:
            records = [
                [handle.readline() for _ in range(4)] for handle in (first, second)
            ]
            if records == [[""] * 4, [""] * 4]:
                return count
            names = []
            for record in records:
                header, sequence, plus, quality = [
                    line.rstrip("\r\n") for line in record
                ]
                if (
                    not header.startswith("@")
                    or not plus.startswith("+")
                    or not sequence
                ):
                    raise ValueError("Incomplete or invalid FASTQ record")
                if len(sequence) != len(quality):
                    raise ValueError("FASTQ sequence and quality lengths differ")
                name = header.split()[0][1:]
                names.append(name[:-2] if name.endswith(("/1", "/2")) else name)
            if names[0] != names[1]:
                raise ValueError("Pilot FASTQ mates have different names")
            count += 1


def validate_trim(current: dict, baseline: dict) -> dict:
    """Verify existing input checksums and complete paired-read accounting."""
    for key in (
        "sample_id",
        "raw_pairs_expected",
        "pairs_entering_fastp",
        "pairs_after_fastp",
    ):
        if current[key] != baseline[key]:
            raise ValueError(
                f"Regenerated trimming differs from the completed baseline: {key}"
            )

    def records(data: dict) -> dict:
        return {r["path"]: (r["size_bytes"], r["sha256"]) for r in data["inputs"]}

    if records(current) != records(baseline):
        raise ValueError(
            "Raw-file checksums differ from the completed fish-screen inputs"
        )
    if current["raw_pairs_expected"] != current["pairs_entering_fastp"]:
        raise ValueError(
            "Halistemma must use all raw pairs; neither library needs the 200M cap"
        )
    return {
        "status": "PASS",
        "sample_id": current["sample_id"],
        "pairs": current["pairs_after_fastp"],
    }


def pilot_metrics(
    sample: str, pairs: int, primary: int, mapped: int, both_unmapped: int
) -> dict:
    if primary != 2 * pairs or both_unmapped % 2 or not 0 <= mapped <= primary:
        raise ValueError("Pilot BAM and FASTQ read counts disagree")
    retained = both_unmapped // 2
    if not 0 <= retained <= pairs or mapped > primary - both_unmapped:
        raise ValueError("Invalid paired mapping counts")
    return {
        "sample_id": sample,
        "input_pairs": pairs,
        "primary_reads": primary,
        "mapped_reads": mapped,
        "mapped_read_fraction": mapped / primary,
        "both_unmapped_pairs": retained,
        "retained_pair_fraction": retained / pairs,
        "host_removed_pair_fraction": 1 - retained / pairs,
    }


def validate_pilots(reports: list[dict], samples: set[str], minimum: float) -> dict:
    if len(reports) != len(samples) or {r["sample_id"] for r in reports} != samples:
        raise ValueError(
            "Pilot reports do not cover exactly the two Halistemma libraries"
        )
    poor = [
        r["sample_id"] for r in reports if r["host_removed_pair_fraction"] < minimum
    ]
    if poor:
        raise ValueError(
            f"Low host mapping requires review before full depletion: {poor}"
        )
    return {
        "status": "PASS",
        "minimum_host_removed_pair_fraction": minimum,
        "libraries": reports,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    trim = sub.add_parser("trim")
    trim.add_argument("--current", required=True, type=Path)
    trim.add_argument("--baseline", required=True, type=Path)
    pilot = sub.add_parser("pilot")
    for key in ("r1", "r2", "bam"):
        pilot.add_argument(f"--{key}", required=True, type=Path)
    pilot.add_argument("--samtools", required=True)
    pilot.add_argument("--sample", required=True)
    pilot.add_argument("--pairs", required=True, type=int)
    gate = sub.add_parser("gate")
    gate.add_argument("--reports", nargs="+", required=True, type=Path)
    gate.add_argument("--samples", required=True, type=Path)
    gate.add_argument("--minimum", required=True, type=float)
    member = sub.add_parser("membership")
    member.add_argument("--eligibility", required=True, type=Path)
    member.add_argument("--samples", required=True, type=Path)
    member.add_argument("--minimum", required=True, type=int)
    for command in (trim, pilot, gate, member):
        command.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "trim":
        result = validate_trim(
            json.loads(args.current.read_text()), json.loads(args.baseline.read_text())
        )
    elif args.command == "pilot":
        count = paired_count(args.r1, args.r2)
        if count != args.pairs:
            raise ValueError(
                f"Pilot must contain exactly {args.pairs} paired reads; found {count}"
            )
        counts = []
        commands = []
        for flags in (["-F", "2304"], ["-F", "2308"], ["-f", "12", "-F", "2304"]):
            command = [args.samtools, "view", "-c", *flags, str(args.bam)]
            counts.append(int(subprocess.check_output(command, text=True).strip()))
            commands.append(command)
        result = pilot_metrics(args.sample, count, *counts)
        result["count_commands"] = commands
        result["samtools_version"] = subprocess.check_output(
            [args.samtools, "--version"], text=True
        ).splitlines()[0]
    else:
        with args.samples.open() as handle:
            samples = {r["sample_id"] for r in csv.DictReader(handle, delimiter="\t")}
        if args.command == "gate":
            result = validate_pilots(
                [json.loads(p.read_text()) for p in args.reports], samples, args.minimum
            )
        else:
            with args.eligibility.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if len(rows) != len(samples) or {r["sample_id"] for r in rows} != samples:
                raise ValueError(
                    "Assembly eligibility does not match the two input libraries"
                )
            if any(
                r["assembly_eligible"] != "true"
                or int(r["retained_pairs"]) < args.minimum
                for r in rows
            ):
                raise ValueError(
                    "Host-depleted input falls below the established assembly floor; review required"
                )
            fields = [
                "sample_id",
                "host_route",
                "strategy",
                "retained_pairs",
                "minimum_assembly_input_pairs",
                "assembly_eligible",
            ]
            args.output.parent.mkdir(parents=True, exist_ok=True)
            with args.output.open("w") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
                )
                writer.writeheader()
                writer.writerows(
                    dict(row, minimum_assembly_input_pairs=args.minimum) for row in rows
                )
            return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
