#!/usr/bin/env python3
"""Compare pilot candidate sets and export aligner-only pairs for verification."""

import argparse
import csv
import json
import shutil
from pathlib import Path

from pilot import digest, pair_ids, read_id, records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--screen-root", type=Path, required=True)
    parser.add_argument("--aligner-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    with args.manifest.open() as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    comparison, selected = [], []
    for row in manifest:
        relative = Path(row["analysis_sample_id"], row["read_pair_id"])
        screen, aligner = args.screen_root / relative, args.aligner_root / relative
        for directory in [screen, aligner]:
            if json.loads((directory / "run.json").read_text())["status"] != "complete":
                raise ValueError(f"Incomplete search: {directory}")
        bb = pair_ids(screen / "bbduk25_R1.fastq", screen / "bbduk25_R2.fastq")
        regex = {
            read_id(r[0]) for m in (1, 2) for r in records(screen / f"regex_R{m}.fastq")
        }
        bowtie = set()
        with (aligner / "bowtie2.sam").open() as handle:
            for line in handle:
                if not line.startswith("@"):
                    fields = line.split("\t")
                    if not int(fields[1]) & 4:
                        bowtie.add(read_id(fields[0]))
        retained = list(records(aligner / "sortmerna.fq"))
        sortmerna = {read_id(r[0]) for r in retained}
        if bowtie - sortmerna:
            raise ValueError("Bowtie2-only pairs need extraction from the input FASTQs")
        new = (bowtie | sortmerna) - (bb | regex)
        comparison.append(
            {
                "analysis_sample_id": row["analysis_sample_id"],
                "bowtie2_equals_sortmerna": bowtie == sortmerna,
                "aligner_only_pairs": len(new),
                "bbduk_only_pairs": len(bb - (bowtie | sortmerna)),
            }
        )
        if not new:
            continue
        selected.append(row)
        output = args.output / relative
        output.mkdir(parents=True)
        counts = {}
        with (
            (output / "candidates.fasta").open("w") as fasta,
            (output / "candidate_reads.tsv").open("w") as table,
        ):
            writer = csv.writer(table, delimiter="\t", lineterminator="\n")
            writer.writerow(["query_id", "read_pair_id", "raw_header"])
            for title, sequence, _ in retained:
                rid = read_id(title)
                if rid in new:
                    counts[rid] = counts.get(rid, 0) + 1
                    query = f"aligner_only_{sum(counts.values()):07d}"
                    fasta.write(f">{query}\n{sequence}\n")
                    writer.writerow([query, rid, title])
        if set(counts) != new or any(count != 2 for count in counts.values()):
            raise ValueError(
                "SortMeRNA output does not contain both mates exactly once"
            )
        shutil.copyfile(__file__, output / Path(__file__).name)
        inputs = [
            screen / "run.json",
            aligner / "run.json",
            aligner / "sortmerna.fq",
            aligner / "bowtie2.sam",
        ]
        (output / "run.json").write_text(
            json.dumps(
                {
                    "status": "complete",
                    "input": row,
                    "new_pairs": len(new),
                    "source_sha256": digest(Path(__file__)),
                    "input_sha256": {str(p): digest(p) for p in inputs},
                },
                indent=2,
            )
            + "\n"
        )
    for name, rows, fields in [
        ("comparison.tsv", comparison, list(comparison[0])),
        ("manifest.tsv", selected, list(manifest[0])),
    ]:
        with (args.output / name).open("w") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    main()
