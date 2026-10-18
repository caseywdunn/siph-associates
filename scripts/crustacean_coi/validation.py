#!/usr/bin/env python3
"""Prepare a held-out-genus sensitivity test and audit preserved pilot candidates."""

import argparse
import csv
import hashlib
import json
import random
import re
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path: Path, rows: list[dict]) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def prepare(config: dict) -> None:
    """Withhold all records mentioning selected genera, including old names."""
    source, root = Path(config["source_panel"]), Path(config["validation_run"])
    panel, reads = root / "panel", root / "heldout"
    panel.mkdir(parents=True, exist_ok=True)
    reads.mkdir(parents=True, exist_ok=True)
    rows = table(source / "reference_metadata.tsv")
    nt = {r.id: str(r.seq) for r in SeqIO.parse(source / "competitive.fna", "fasta")}
    eligible = defaultdict(dict)
    for row in rows:
        words = row["species"].split()
        if (
            row["source"] == "RefSeq"
            and row["target"] == "1"
            and int(row["length"]) >= 1200
            and len(words) == 2
            and re.fullmatch(r"[a-z]+", words[1])
        ):
            eligible[row["group"]].setdefault(words[0], row)
    rng = random.Random(20261017)
    selected = []
    for group, genera in sorted(eligible.items()):
        for genus in rng.sample(sorted(genera), min(5, len(genera))):
            selected.append(genera[genus])
    withheld = {r["species"].split()[0] for r in selected}
    held_sequences = {nt[r["reference_id"]] for r in selected}
    held_sequences |= {str(Seq(s).reverse_complement()) for s in held_sequences}
    excluded = {
        r["reference_id"]
        for r in rows
        if set(re.findall(r"[A-Za-z]+", r["species"])) & withheld
        or nt[r["reference_id"]] in held_sequences
    }
    retained = [r for r in rows if r["reference_id"] not in excluded]
    assert excluded and retained
    write_table(panel / "reference_metadata.tsv", retained)
    for name in ("target.fna", "target.faa", "competitive.fna"):
        SeqIO.write(
            (r for r in SeqIO.parse(source / name, "fasta") if r.id not in excluded),
            panel / name,
            "fasta",
        )
    truth = []
    with (reads / "R1.fastq").open("w") as r1, (reads / "R2.fastq").open("w") as r2:
        for row in selected:
            seq = nt[row["reference_id"]]
            for i in range(5):
                start = (len(seq) - 300) * i // 4
                rid = f"heldout_{len(truth) + 1:04d}"
                fragment = seq[start : start + 300]
                for mate, handle, sequence in (
                    (1, r1, fragment[:150]),
                    (2, r2, str(Seq(fragment[-150:]).reverse_complement())),
                ):
                    handle.write(
                        f"@{rid} {mate}:N:0:validation\n{sequence}\n+\n{'I' * 150}\n"
                    )
                truth.append(
                    dict(
                        read_pair_id=rid,
                        expected="crustacean",
                        group=row["group"],
                        species=row["species"],
                        reference_id=row["reference_id"],
                        start=start,
                    )
                )
        # Host, fish, and insect negatives use the original paired controls.
        controls = Path(config["pilot_run"]) / "controls"
        original = json.loads((controls / "truth.json").read_text())
        negative = {
            r["read_pair_id"] for r in original if r["expected"] != "crustacean"
        }
        for mate, handle in ((1, r1), (2, r2)):
            for record in SeqIO.parse(controls / f"R{mate}.fastq", "fastq"):
                if record.id in negative:
                    SeqIO.write(record, handle, "fastq")
        truth.extend(
            dict(r, group="negative", species=r["read_pair_id"], start="")
            for r in original
            if r["read_pair_id"] in negative
        )
    (reads / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")
    (root / "design.json").write_text(
        json.dumps(
            dict(
                seed=20261017,
                withheld_genera=sorted(withheld),
                source_records=len(rows),
                excluded_records=len(excluded),
                retained_records=len(retained),
                selected_references=selected,
                positive_pairs=len(truth) - len(negative),
                negative_pairs=len(negative),
                limitations="Genus names are source labels, not reconciled taxonomy; five error-free 300-nt fragments per reference, paired 150-nt reads. This tests reference divergence, not abundance or sequencing-error sensitivity.",
            ),
            indent=2,
        )
        + "\n"
    )


def audit_queries(config: dict) -> None:
    root = Path(config["validation_run"]) / "audit"
    root.mkdir(parents=True, exist_ok=True)
    mapping = []
    with (root / "queries.fna").open("w") as handle:
        for directory in sorted(Path(config["pilot_run"]).glob("Ahuja2026__*")):
            pairs = {
                r["candidate_id"]: r for r in table(directory / "candidate_pairs.tsv")
            }
            for record in SeqIO.parse(directory / "candidates.fna", "fasta"):
                cid = record.id.rsplit("_", 1)[0]
                query = f"{directory.name}__{record.id}"
                handle.write(f">{query}\n{record.seq}\n")
                mapping.append(dict(query=query, sample=directory.name, **pairs[cid]))
        # Include the barcodes responsible for the initial NA35 labels.
        suspect = {
            "MZG_Amphipoda_2488",
            "MZG_Amphipoda_6928",
            "MZG_Amphipoda_6933",
            "MZG_Decapoda_56382",
        }
        for record in SeqIO.parse(Path(config["source_panel"]) / "target.fna", "fasta"):
            if record.id in suspect:
                handle.write(f">REF__{record.id}\n{record.seq}\n")
    write_table(root / "query_pairs.tsv", mapping)


def report(config: dict) -> None:
    root = Path(config["validation_run"])
    truth = json.loads((root / "heldout/truth.json").read_text())
    observed = {
        r["read_pair_id"]: r for r in table(root / "heldout/classified_pairs.tsv")
    }
    counts = defaultdict(Counter)
    details = []
    for row in truth:
        obs = observed.get(row["read_pair_id"], {})
        detail = dict(
            row,
            classification=obs.get("classification", "not_recruited"),
            bbduk=int(obs.get("bbduk", 0)),
            diamond=int(obs.get("diamond", 0)),
            assigned_groups=obs.get("groups", ""),
        )
        details.append(detail)
        count = counts[row["group"]]
        count["pairs"] += 1
        count["bbduk"] += detail["bbduk"]
        count["diamond"] += detail["diamond"]
        count["union"] += bool(obs)
        count[detail["classification"]] += 1
    write_table(root / "sensitivity_pairs.tsv", details)
    (root / "sensitivity_summary.json").write_text(json.dumps(counts, indent=2) + "\n")
    # Retain the leading broad-reference matches, including all exact score ties.
    for filename, score_column in (
        ("mitochondrial.tsv", 11),
        ("swissprot_code1.tsv", 11),
        ("swissprot_code5.tsv", 11),
    ):
        leading = defaultdict(list)
        with (root / "audit" / filename).open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                query, score = fields[0], float(fields[score_column])
                if not leading[query] or score > float(leading[query][0][score_column]):
                    leading[query] = [fields]
                elif score == float(leading[query][0][score_column]):
                    leading[query].append(fields)
        with (root / "audit" / filename.replace(".tsv", "_best.tsv")).open(
            "w"
        ) as handle:
            for hits in leading.values():
                for fields in hits:
                    handle.write("\t".join(fields) + "\n")
    with (root / "output_sha256.txt").open("w") as handle:
        for path in sorted(root.rglob("*.tsv")) + sorted(root.rglob("*.json")):
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["prepare", "audit", "report"])
    parser.add_argument(
        "--config", type=Path, default=Path("config/crustacean_coi_validation.json")
    )
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    {"prepare": prepare, "audit": audit_queries, "report": report}[args.stage](config)


if __name__ == "__main__":
    main()
