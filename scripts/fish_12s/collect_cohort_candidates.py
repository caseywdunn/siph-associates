#!/usr/bin/env python3
"""Retain both original mates for the union of the four tested fish filters."""

import argparse
import csv
import json
from itertools import zip_longest
from pathlib import Path

from pilot import digest, read_id, records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["manifest", "screen", "aligners", "output"]:
        parser.add_argument(f"--{name}", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    with args.manifest.open() as handle:
        row = next(
            r
            for r in csv.DictReader(handle, delimiter="\t")
            if r["analysis_sample_id"] == args.sample
        )
    for path in [args.screen / "run.json", args.aligners / "run.json"]:
        if json.loads(path.read_text())["status"] != "complete":
            raise ValueError(f"Incomplete search: {path}")
    selected = {
        "bbduk25": {read_id(r[0]) for r in records(args.screen / "bbduk25_R1.fastq")},
        "regex": {
            read_id(r[0])
            for mate in [1, 2]
            for r in records(args.screen / f"regex_R{mate}.fastq")
        },
        "sortmerna": {read_id(r[0]) for r in records(args.aligners / "sortmerna.fq")},
        "bowtie2": set(),
    }
    with (args.aligners / "bowtie2.sam").open() as handle:
        for line in handle:
            if line.startswith("@"):
                continue
            fields = line.split("\t")
            if not int(fields[1]) & 4:
                selected["bowtie2"].add(read_id(fields[0]))
    wanted, seen = set.union(*selected.values()), set()
    args.output.mkdir(parents=True, exist_ok=True)
    count = 0
    with (
        (args.output / "candidates.fasta").open("w") as fasta,
        (args.output / "candidate_pairs.tsv").open("w") as table,
    ):
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(["candidate_id", "read_pair_id", *selected])
        for first, second in zip_longest(
            records(Path(row["r1_path"])), records(Path(row["r2_path"]))
        ):
            if (
                first is None
                or second is None
                or read_id(first[0]) != read_id(second[0])
            ):
                raise ValueError("Unsynchronized retained reads")
            count += 1
            rid = read_id(first[0])
            if rid not in wanted:
                continue
            if rid in seen:
                raise ValueError(f"Duplicate candidate pair: {rid}")
            seen.add(rid)
            cid = f"c{len(seen):08d}"
            writer.writerow([cid, rid, *(int(rid in ids) for ids in selected.values())])
            for mate, record in enumerate([first, second], 1):
                fasta.write(f">{cid}_{mate}\n{record[1]}\n")
    if count != int(row["pairs"]) or seen != wanted:
        raise ValueError("Fish candidate/input pair accounting failed")
    (args.output / "run.json").write_text(
        json.dumps(
            dict(
                status="complete",
                input=row,
                input_pairs=count,
                union_pairs=len(seen),
                method_pairs={name: len(ids) for name, ids in selected.items()},
                pair_synchronization_checked=True,
                search_records={
                    str(path): digest(path)
                    for path in [args.screen / "run.json", args.aligners / "run.json"]
                },
                source_sha256={
                    p.name: digest(p)
                    for p in [Path(__file__), Path(__file__).with_name("pilot.py")]
                },
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
