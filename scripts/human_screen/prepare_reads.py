"""Select per-source-file FASTQ prefixes or simulate known-origin control pairs."""

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import random
from pathlib import Path

import pysam


def read_record(handle):
    """Read a validated four-line FASTQ record, preserving its text."""
    header = handle.readline()
    if not header:
        return None
    seq, plus, qual = (handle.readline().rstrip("\n") for _ in range(3))
    if not header.startswith("@") or not plus.startswith("+") or len(seq) != len(qual):
        raise ValueError(f"Malformed FASTQ record: {header.strip()}")
    return header, seq, plus, qual


def read_name(header):
    return header.split()[0].removesuffix("/1").removesuffix("/2")


def write_record(handle, name, sequence, quality):
    handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


def prepare_prefix(row, pairs, out1, out2):
    """Take equal-sized prefixes from all original source files."""
    paths1, paths2 = row["r1_paths"].split(";"), row["r2_paths"].split(";")
    if len(paths1) != len(paths2) or len(set(paths1)) != len(paths1):
        raise ValueError("Mismatched or repeated source files")
    inventory = []
    with (
        gzip.open(out1, "wt", compresslevel=1) as w1,
        gzip.open(out2, "wt", compresslevel=1) as w2,
    ):
        for lane, (path1, path2) in enumerate(zip(paths1, paths2, strict=True)):
            target = pairs // len(paths1) + int(lane < pairs % len(paths1))
            digest = hashlib.sha256()
            count = 0
            with gzip.open(path1, "rt") as h1, gzip.open(path2, "rt") as h2:
                for index in range(target):
                    r1, r2 = read_record(h1), read_record(h2)
                    if r1 is None or r2 is None:
                        if r1 != r2:
                            raise ValueError("Unequal mate counts")
                        break
                    if read_name(r1[0]) != read_name(r2[0]):
                        raise ValueError(
                            f"Unpaired FASTQ identifiers: {r1[0]} / {r2[0]}"
                        )
                    # Source index prevents repeated lane names from colliding downstream.
                    name = f"L{lane}_{index}"
                    write_record(w1, name, r1[1], r1[3])
                    write_record(w2, name, r2[1], r2[3])
                    digest.update((r1[1] + r2[1]).encode())
                    count += 1
            inventory.append(
                {
                    "r1": path1,
                    "r2": path2,
                    "requested_pairs": target,
                    "selected_pairs": count,
                    "selected_sequence_sha256": digest.hexdigest(),
                    "r1_bytes": Path(path1).stat().st_size,
                    "r2_bytes": Path(path2).stat().st_size,
                }
            )
    return {
        "selection": "equal_per_source_file_prefix",
        "source_files": inventory,
        "selected_pairs": sum(x["selected_pairs"] for x in inventory),
    }


def simulate(reference, count, seed, out1, out2):
    """Simulate 150-bp reads from 350-bp inserts, one substitution per mate."""
    rng = random.Random(seed)
    # Build a local FASTA index if the shared reference lacks one; never mutate it.
    index = str(out1) + ".reference.fai"
    if Path(str(reference) + ".fai").exists():
        fasta = pysam.FastaFile(str(reference))
    else:
        pysam.faidx(str(reference), "--fai-idx", index)
        fasta = pysam.FastaFile(str(reference), filepath_index=index)
    eligible = [
        (name, length - 350 + 1)
        for name, length in zip(fasta.references, fasta.lengths, strict=True)
        if length >= 350
    ]
    cumulative, total = [], 0
    for _, length in eligible:
        total += length
        cumulative.append(total)
    comp = str.maketrans("ACGT", "TGCA")
    emitted = 0
    with (
        gzip.open(out1, "wt", compresslevel=1) as w1,
        gzip.open(out2, "wt", compresslevel=1) as w2,
    ):
        while emitted < count:
            draw = rng.randrange(total)
            chrom_index = bisect.bisect_right(cumulative, draw)
            offset = draw - (cumulative[chrom_index - 1] if chrom_index else 0)
            fragment = fasta.fetch(
                eligible[chrom_index][0], offset, offset + 350
            ).upper()
            if set(fragment) - set("ACGT"):
                continue
            sequences = [fragment[:150], fragment[-150:].translate(comp)[::-1]]
            for mate in range(2):
                seq = list(sequences[mate])
                pos = rng.randrange(150)
                seq[pos] = rng.choice([b for b in "ACGT" if b != seq[pos]])
                sequences[mate] = "".join(seq)
            write_record(w1, f"control_{emitted}", sequences[0], "I" * 150)
            write_record(w2, f"control_{emitted}", sequences[1], "I" * 150)
            emitted += 1
    fasta.close()
    return {
        "selection": "simulated_reference_pairs",
        "reference": str(reference),
        "seed": seed,
        "selected_pairs": emitted,
        "read_length": 150,
        "insert_length": 350,
        "substitutions_per_mate": 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", required=True, type=Path)
    parser.add_argument("--r2", required=True, type=Path)
    parser.add_argument("--metrics", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text())
    args.r1.parent.mkdir(parents=True, exist_ok=True)
    if args.sample.startswith("control_"):
        label = args.sample.removeprefix("control_")
        reference = (
            config["human_fasta"]
            if label == "human"
            else config["negative_references"][label]
        )
        metrics = simulate(
            reference, config["control_pairs"], config["seed"], args.r1, args.r2
        )
    else:
        with open(config["manifest"]) as handle:
            rows = [
                row
                for row in csv.DictReader(handle, delimiter="\t")
                if row["sample_id"] == args.sample
            ]
        if len(rows) != 1:
            raise ValueError(f"Expected one manifest entry for {args.sample}")
        metrics = prepare_prefix(rows[0], config["pairs_per_library"], args.r1, args.r2)
    metrics["sample_id"] = args.sample
    args.metrics.write_text(json.dumps(metrics, indent=2) + "\n")


if __name__ == "__main__":
    main()
