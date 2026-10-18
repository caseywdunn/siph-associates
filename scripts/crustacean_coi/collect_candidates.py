#!/usr/bin/env python3
"""Keep both mates for the union of nucleotide and translated COI matches."""

import argparse
import csv
import gzip
import json
from itertools import zip_longest
from pathlib import Path

from Bio.SeqIO.QualityIO import FastqGeneralIterator


def pair_id(title):
    value = title.split()[0]
    return value[:-2] if value.endswith(("/1", "/2")) else value


def fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        yield from FastqGeneralIterator(handle)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--r1", type=Path, required=True)
    parser.add_argument("--r2", type=Path, required=True)
    parser.add_argument("--bbduk", type=Path, required=True)
    parser.add_argument("--diamond", type=Path, nargs=2, required=True)
    parser.add_argument("--expected-pairs", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    nucleotide = {pair_id(r[0]) for r in fastq(args.bbduk)}
    translated = set()
    for path in args.diamond:
        with path.open() as handle:
            translated.update(
                pair_id(line.split()[0]) for line in handle if line.strip()
            )
    wanted = nucleotide | translated
    seen, count = set(), 0
    args.output.mkdir(parents=True, exist_ok=True)
    with (
        (args.output / "candidates.fna").open("w") as fasta,
        (args.output / "candidate_pairs.tsv").open("w") as table,
    ):
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(["candidate_id", "read_pair_id", "bbduk", "diamond"])
        for r1, r2 in zip_longest(fastq(args.r1), fastq(args.r2)):
            if r1 is None or r2 is None or pair_id(r1[0]) != pair_id(r2[0]):
                raise ValueError(
                    f"Input mates are not synchronized at pair {count + 1}"
                )
            count += 1
            rid = pair_id(r1[0])
            if rid not in wanted:
                continue
            if rid in seen:
                raise ValueError(f"Duplicate candidate identifier: {rid}")
            seen.add(rid)
            cid = f"c{len(seen):08d}"
            writer.writerow([cid, rid, int(rid in nucleotide), int(rid in translated)])
            for mate, record in enumerate([r1, r2], 1):
                fasta.write(f">{cid}_{mate}\n{record[1]}\n")
    if count != args.expected_pairs or seen != wanted:
        raise ValueError(
            f"Pair accounting failed: {count}/{args.expected_pairs}; {len(wanted - seen)} missing candidates"
        )
    (args.output / "recruitment.json").write_text(
        json.dumps(
            {
                "input_pairs": count,
                "bbduk_pairs": len(nucleotide),
                "diamond_pairs": len(translated),
                "both_methods_pairs": len(nucleotide & translated),
                "union_pairs": len(wanted),
                "pair_synchronization_checked": True,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
