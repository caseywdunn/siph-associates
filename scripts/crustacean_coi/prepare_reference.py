#!/usr/bin/env python3
"""Extract annotated COI and validate the reading frames of public barcodes."""

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq, UndefinedSequenceError

GROUPS = (
    "Copepoda",
    "Euphausiacea",
    "Decapoda",
    "Amphipoda",
    "Mysida",
    "Lophogastrida",
    "Ostracoda",
)
FIELDS = [
    "reference_id",
    "source",
    "source_id",
    "species",
    "group",
    "target",
    "genetic_code",
    "length",
    "source_description",
]


def is_coi(feature) -> bool:
    """Recognize COI without accidentally accepting COII or COIII."""
    if feature.type != "CDS":
        return False
    genes = {
        re.sub(r"[^a-z0-9]", "", x.lower()) for x in feature.qualifiers.get("gene", [])
    }
    if genes & {"coi", "co1", "cox1", "coxi", "mtco1"}:
        return True
    return any(
        re.search(r"cytochrome\s+c\s+oxidase\s+(?:subunit\s+)?(?:i|1)\b", x, re.I)
        for x in feature.qualifiers.get("product", [])
    )


def translate_frame(sequence: str, frame: int, code: int = 5) -> str:
    if frame not in (-3, -2, -1, 1, 2, 3):
        raise ValueError(f"Invalid reading frame: {frame}")
    oriented = Seq(sequence) if frame > 0 else Seq(sequence).reverse_complement()
    coding = oriented[abs(frame) - 1 :]
    return str(coding[: len(coding) // 3 * 3].translate(table=code)).rstrip("*")


def table(path: Path, rows: list[dict]) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def extract(args) -> None:
    args.output.mkdir(parents=True, exist_ok=True)
    rows, excluded = [], Counter()
    sources = json.loads(args.sources.read_text())
    with (
        (args.output / "annotated.fna").open("w") as nt,
        (args.output / "annotated.faa").open("w") as aa,
    ):
        source = next(s for s in sources if s["kind"] == "genbank")
        with gzip.open(args.reference_dir / source["filename"], "rt") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                lineage = record.annotations.get("taxonomy", [])
                group = next((g for g in GROUPS if g in lineage), "non_target")
                if "Cnidaria" in lineage:
                    group = "Cnidaria"
                for index, feature in enumerate(record.features):
                    if not is_coi(feature):
                        continue
                    try:
                        sequence = str(feature.extract(record.seq)).upper()
                    except UndefinedSequenceError:
                        # Some RefSeq entries describe a contig without supplying
                        # its bases in this release file.
                        excluded["annotated_undefined_sequence"] += 1
                        continue
                    protein = feature.qualifiers.get("translation", [""])[0]
                    code = int(feature.qualifiers.get("transl_table", ["1"])[0])
                    offset = int(feature.qualifiers.get("codon_start", ["1"])[0])
                    if not 180 <= len(sequence) <= 2400 or set(sequence) - set("ACGT"):
                        excluded["annotated_length_or_ambiguity"] += 1
                        continue
                    translated = translate_frame(sequence, offset, code)
                    # Alternative initiation codons are translated as M in
                    # annotated CDSs, but as their ordinary amino acid here.
                    valid = (
                        protein and "*" not in protein and translated[1:] == protein[1:]
                    )
                    if group in GROUPS and (not valid or code != 5):
                        excluded[f"target_translation_or_code_{code}"] += 1
                        continue
                    rid = f"GB_{record.id}_{index}"
                    nt.write(f">{rid}\n{sequence}\n")
                    if valid:
                        aa.write(f">{rid}\n{protein}\n")
                    rows.append(
                        dict(
                            zip(
                                FIELDS,
                                [
                                    rid,
                                    "RefSeq",
                                    record.id,
                                    record.annotations.get(
                                        "organism", record.description
                                    ),
                                    group,
                                    int(group in GROUPS),
                                    code,
                                    len(sequence),
                                    record.description,
                                ],
                            )
                        )
                    )
    table(args.output / "annotated_metadata.tsv", rows)
    barcodes = []
    with (
        (args.output / "barcodes.fna").open("w") as nt,
        (args.output / "host_barcodes.fna").open("w") as host,
    ):
        for source in (s for s in sources if s["kind"] == "mzg_coi"):
            with gzip.open(args.reference_dir / source["filename"], "rt") as handle:
                for index, record in enumerate(SeqIO.parse(handle, "fasta")):
                    sequence = str(record.seq).upper()
                    if not 180 <= len(sequence) <= 2400 or set(sequence) - set("ACGT"):
                        excluded[f"{source['name']}_length_or_ambiguity"] += 1
                        continue
                    rid = f"MZG_{source['name']}_{index}"
                    species = record.id.partition("__")[2].replace("_", " ")
                    if not species:
                        raise ValueError(
                            f"Unrecognized MetaZooGene header: {record.description}"
                        )
                    (nt if source["target"] else host).write(f">{rid}\n{sequence}\n")
                    barcodes.append(
                        dict(
                            zip(
                                FIELDS,
                                [
                                    rid,
                                    "MetaZooGene",
                                    record.id,
                                    species,
                                    source["name"],
                                    int(source["target"]),
                                    5 if source["target"] else "",
                                    len(sequence),
                                    record.description,
                                ],
                            )
                        )
                    )
    table(args.output / "barcode_metadata.tsv", barcodes)
    (args.output / "extraction.json").write_text(
        json.dumps(
            {
                "annotated_records": len(rows),
                "barcodes": len(barcodes),
                "excluded": excluded,
                "annotated_target_groups": dict(
                    Counter(r["group"] for r in rows if r["target"])
                ),
            },
            indent=2,
        )
        + "\n"
    )


def finalize(args) -> None:
    out = args.output
    annotated = list(
        csv.DictReader((out / "annotated_metadata.tsv").open(), delimiter="\t")
    )
    barcodes = list(
        csv.DictReader((out / "barcode_metadata.tsv").open(), delimiter="\t")
    )
    rows = {r["reference_id"]: r for r in annotated + barcodes}
    nt = {
        r.id: str(r.seq)
        for name in ["annotated.fna", "barcodes.fna", "host_barcodes.fna"]
        for r in SeqIO.parse(out / name, "fasta")
    }
    aa = {r.id: str(r.seq) for r in SeqIO.parse(out / "annotated.faa", "fasta")}
    frames = {}
    with (out / "barcode_frames.tsv").open() as handle:
        for line in handle:
            rid, _, identity, length, coverage, frame, score = line.split()
            if (
                float(identity) >= 50
                and int(length) >= 50
                and float(coverage) >= 80
                and (rid not in frames or float(score) > frames[rid][0])
            ):
                frames[rid] = (float(score), int(frame))
    rejected = Counter()
    accepted = []
    for row in barcodes:
        rid = row["reference_id"]
        if row["target"] != "1":
            continue
        if rid not in frames:
            rejected["no_supported_coi_frame"] += 1
            continue
        protein = translate_frame(nt[rid], frames[rid][1])
        if "*" in protein or len(protein) < 50:
            rejected["stop_or_short_translation"] += 1
            continue
        aa[rid] = protein
        accepted.append(rid)
    targets = {r["reference_id"] for r in annotated if r["target"] == "1"} | set(
        accepted
    )
    competitors = {
        r["reference_id"] for r in annotated + barcodes if r["target"] != "1"
    }
    if not targets or any(rid not in aa for rid in targets):
        raise ValueError("Target nucleotide and protein panels are inconsistent")
    with (
        (out / "target.fna").open("w") as fna,
        (out / "target.faa").open("w") as faa,
        (out / "competitive.fna").open("w") as comp,
    ):
        for rid in sorted(targets | competitors):
            comp.write(f">{rid}\n{nt[rid]}\n")
            if rid in targets:
                fna.write(f">{rid}\n{nt[rid]}\n")
                faa.write(f">{rid}\n{aa[rid]}\n")
    table(
        out / "reference_metadata.tsv",
        [rows[rid] for rid in sorted(targets | competitors)],
    )
    coverage = []
    for group in GROUPS:
        selected = [rows[r] for r in targets if rows[r]["group"] == group]
        coverage.append(
            {
                "group": group,
                "sequences": len(selected),
                "species_labels": len({r["species"] for r in selected}),
                "RefSeq": sum(r["source"] == "RefSeq" for r in selected),
                "MetaZooGene": sum(r["source"] == "MetaZooGene" for r in selected),
            }
        )
    with (out / "coverage.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle, list(coverage[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(coverage)
    (out / "reference_summary.json").write_text(
        json.dumps(
            {
                "target_records": len(targets),
                "competitor_records": len(competitors),
                "rejected_barcodes": rejected,
                "frame_validation": "DIAMOND >=50aa, >=50% identity, >=80% query coverage; no internal stop; code 5",
                "barcode_frames": {r: frames[r][1] for r in accepted},
            },
            indent=2,
        )
        + "\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["extract", "finalize"])
    parser.add_argument("--reference-dir", type=Path)
    parser.add_argument("--sources", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    (extract if args.stage == "extract" else finalize)(args)


if __name__ == "__main__":
    main()
