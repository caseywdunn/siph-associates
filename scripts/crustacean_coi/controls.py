#!/usr/bin/env python3
"""Generate and check paired COI controls, including conserved host sequences."""

import argparse
import csv
import json
import random
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq


def generate(panel, output):
    output.mkdir(parents=True, exist_ok=True)
    with (panel / "reference_metadata.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sequences = {
        r.id: str(r.seq) for r in SeqIO.parse(panel / "competitive.fna", "fasta")
    }
    rng = random.Random(20261017)

    def random_read():
        return "".join(rng.choices("ACGT", k=150))

    records = []
    for group in ["Copepoda", "Euphausiacea", "Decapoda"]:
        row = next(r for r in rows if r["group"] == group and r["source"] == "RefSeq")
        seq = sequences[row["reference_id"]][300:450]
        records.append(
            (f"exact_{group}", seq, random_read(), "crustacean", row["reference_id"])
        )
        if group == "Copepoda":
            records.append(
                (
                    "reverse_R2_copepod",
                    random_read(),
                    str(Seq(seq).reverse_complement()),
                    "crustacean",
                    row["reference_id"],
                )
            )
            code = CodonTable.unambiguous_dna_by_id[5].forward_table
            altered = list(seq)
            for i in range(0, len(seq) - 2, 9):
                codon = seq[i : i + 3]
                synonyms = [
                    c for c, aa in code.items() if aa == code.get(codon) and c != codon
                ]
                if synonyms:
                    altered[i : i + 3] = sorted(synonyms)[0]
            changed = "".join(altered)
            assert Seq(seq).translate(table=5) == Seq(changed).translate(table=5)
            records.append(
                (
                    "synonymous_copepod",
                    changed,
                    random_read(),
                    "crustacean",
                    row["reference_id"],
                )
            )
    for name in [
        "Nanomia",
        "Praya",
        "Physalia",
        "Danio rerio",
        "Drosophila melanogaster",
    ]:
        matching = [
            r for r in rows if r["species"].startswith(name) and r["target"] == "0"
        ]
        if not matching:
            raise ValueError(f"Missing requested non-target control: {name}")
        row = max(matching, key=lambda r: int(r["length"]))
        sequence = sequences[row["reference_id"]]
        start = min(300, len(sequence) - 150)
        records.append(
            (
                name.replace(" ", "_"),
                sequence[start : start + 150],
                random_read(),
                "non_target",
                row["reference_id"],
            )
        )
    records.extend(
        [
            ("random", random_read(), random_read(), "negative", ""),
            ("polyG", "G" * 150, random_read(), "negative", ""),
        ]
    )
    with (output / "R1.fastq").open("w") as r1, (output / "R2.fastq").open("w") as r2:
        for name, seq1, seq2, _, _ in records:
            for mate, handle, seq in [(1, r1, seq1), (2, r2, seq2)]:
                handle.write(
                    f"@{name} {mate}:N:0:control\n{seq}\n+\n{'I' * len(seq)}\n"
                )
    truth = [
        {"read_pair_id": r[0], "expected": r[3], "reference_id": r[4]} for r in records
    ]
    (output / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")


def check(output):
    truth = json.loads((output / "truth.json").read_text())
    with (output / "classified_pairs.tsv").open() as handle:
        observed = {
            r["read_pair_id"]: r for r in csv.DictReader(handle, delimiter="\t")
        }
    results = []
    for row in truth:
        actual = observed.get(row["read_pair_id"], {})
        label = actual.get("classification", "not_recruited")
        if row["expected"] == "crustacean":
            passed = label == "crustacean" and actual.get("diamond") == "1"
            if row["read_pair_id"].startswith(("exact_", "reverse_")):
                passed = passed and actual.get("bbduk") == "1"
        else:
            passed = label in ("not_recruited", "non_target")
        results.append(
            dict(
                row,
                observed=label,
                passed=passed,
                bbduk=actual.get("bbduk", "0"),
                diamond=actual.get("diamond", "0"),
            )
        )
    (output / "control_results.json").write_text(json.dumps(results, indent=2) + "\n")
    if not all(r["passed"] for r in results):
        raise ValueError(
            f"COI controls failed: {[r for r in results if not r['passed']]}"
        )
    (output / "PASSED").write_text(f"{len(results)} paired controls passed\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["generate", "check"])
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.stage == "generate":
        generate(args.panel, args.output)
    else:
        check(args.output)


if __name__ == "__main__":
    main()
