#!/usr/bin/env python3
"""Validate and stage the existing bac120 alignment for an IQ-TREE rerun."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from Bio import Phylo, SeqIO


FOCAL_MAGS = {
    "MAGSP0005", "MAGSP0007", "MAGSP0010", "MAGSP0011",
    "MAGSP0012", "MAGSP0029", "MAGSP0031",
}
CANONICAL_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
VALID_AMINO_ACIDS = CANONICAL_AMINO_ACIDS | set("BXZJUO?-")
OUTGROUP_TAXON = "f__Acholeplasmataceae"


def sha256(path: Path) -> str:
    """Return the checksum of a source or staged file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    """Write a compact table with stable column order."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def prepare(alignment: Path, tree_path: Path, labels: Path, output: Path) -> dict:
    """Validate this fixed analysis input before writing any derivatives."""
    opener = gzip.open if alignment.suffix == ".gz" else open
    with opener(alignment, "rt") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    identifiers = [record.id for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("Alignment identifiers are not unique")
    if len(records) != 483:
        raise ValueError(f"Expected 483 sequences; found {len(records)}")
    lengths = {len(record.seq) for record in records}
    if lengths != {5035}:
        raise ValueError(f"Expected 5,035 columns in every sequence; found {lengths}")
    identifier_set = set(identifiers)
    if identifier_set.intersection(FOCAL_MAGS) != FOCAL_MAGS:
        raise ValueError(f"Missing focal genomes: {sorted(FOCAL_MAGS - identifier_set)}")
    if {name for name in identifiers if name.startswith("MAGSP")} != FOCAL_MAGS:
        raise ValueError("Focal MAG set differs from the seven original representatives")

    sequence_rows = []
    sequence_groups: dict[str, list[str]] = defaultdict(list)
    accession_groups: dict[str, list[str]] = defaultdict(list)
    for record in records:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", record.id):
            raise ValueError(f"Unsupported IQ-TREE identifier: {record.id!r}")
        sequence = str(record.seq)
        invalid = set(sequence) - VALID_AMINO_ACIDS
        if invalid:
            raise ValueError(f"Invalid amino-acid characters in {record.id}: {invalid}")
        resolved = sum(character in CANONICAL_AMINO_ACIDS for character in sequence)
        if not resolved:
            raise ValueError(f"No resolved amino-acid positions in {record.id}")
        gaps = sequence.count("-")
        unknown = sum(sequence.count(character) for character in "X?")
        sequence_digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        sequence_groups[sequence_digest].append(record.id)
        accession_match = re.search(r"GC[AF]_\d+\.\d+", record.id)
        accession = accession_match.group() if accession_match else ""
        if accession:
            accession_groups[accession].append(record.id)
        sequence_rows.append({
            "genome": record.id,
            "alignment_columns": len(sequence),
            "resolved_amino_acids": resolved,
            "gap_columns": gaps,
            "unknown_columns": unknown,
            "other_ambiguous_columns": len(sequence) - resolved - gaps - unknown,
            "missing_fraction": (gaps + unknown) / len(sequence),
            "assembly_accession": accession,
            "aligned_sequence_sha256": sequence_digest,
        })

    tree = Phylo.read(tree_path, "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    if len(tips) != len(set(tips)) or set(tips) != identifier_set:
        raise ValueError("Tree tips must be unique and match the alignment exactly")
    outgroup_children = [
        child for child in tree.root.clades
        if OUTGROUP_TAXON in re.split(r"[:;]\s*", child.name or "")
    ]
    if len(tree.root.clades) != 2 or len(outgroup_children) != 1:
        raise ValueError("The original tree does not identify one root-child outgroup")
    outgroup_clade = outgroup_children[0]
    outgroup = sorted(tip.name for tip in outgroup_clade.get_terminals())
    if not outgroup or len(outgroup) == len(tips) or FOCAL_MAGS.intersection(outgroup):
        raise ValueError("Invalid Acholeplasmataceae outgroup in the original tree")

    with labels.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome", "group", "classification", "status"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing label-table columns: {required}")
        label_rows = list(reader)
    label_ids = [row["genome"] for row in label_rows]
    if len(label_ids) != len(set(label_ids)):
        raise ValueError("Duplicate genome identifiers in the label table")
    included = {row["genome"] for row in label_rows if row["status"] == "included"}
    if not included.issubset(identifier_set) or not FOCAL_MAGS.issubset(included):
        raise ValueError("Included label-table genomes do not match the alignment")

    duplicates = []
    for kind, groups in (
        ("identical_aligned_sequence", sequence_groups),
        ("same_assembly_accession", accession_groups),
    ):
        for key, members in sorted(groups.items()):
            if len(members) > 1:
                for member in sorted(members):
                    duplicates.append({
                        "duplicate_kind": kind, "group_key": key,
                        "group_size": len(members), "genome": member,
                    })

    output.mkdir(parents=True, exist_ok=True)
    # Retain every tip; duplicate metadata describe the original analysis input.
    SeqIO.write(records, output / "alignment.fasta", "fasta")
    write_tsv(output / "sequence_summary.tsv", list(sequence_rows[0]), sequence_rows)
    write_tsv(output / "duplicate_sequences.tsv",
              ["duplicate_kind", "group_key", "group_size", "genome"], duplicates)
    (output / "outgroup.txt").write_text(",".join(outgroup) + "\n")
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "sources": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in (("alignment", alignment), ("tree", tree_path), ("labels", labels))
        },
        "sequence_count": len(records),
        "alignment_columns": 5035,
        "focal_genomes": sorted(FOCAL_MAGS),
        "included_label_count": len(included),
        "outgroup_taxon": OUTGROUP_TAXON,
        "outgroup_root_child_label": outgroup_clade.name,
        "outgroup_tip_count": len(outgroup),
        "outgroup_tips": outgroup,
        "identical_sequence_groups": sum(len(group) > 1 for group in sequence_groups.values()),
        "duplicate_accession_groups": sum(len(group) > 1 for group in accession_groups.values()),
        "duplicate_policy": "retain all original alignment sequences and tree tips",
        "outputs": {
            name: sha256(output / name)
            for name in ("alignment.fasta", "sequence_summary.tsv", "duplicate_sequences.tsv", "outgroup.txt")
        },
    }
    (output / "source.sha256").write_text("".join(
        f"{source['sha256']}  {source['path']}\n"
        for source in report["sources"].values()
    ))
    (output / "generated.sha256").write_text("".join(
        f"{digest}  {name}\n" for name, digest in report["outputs"].items()
    ))
    (output / "inputs.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.alignment, args.tree, args.labels, args.output_dir)
    print(json.dumps({key: report[key] for key in (
        "sequence_count", "alignment_columns", "outgroup_tip_count",
        "identical_sequence_groups", "duplicate_accession_groups",
    )}))


if __name__ == "__main__":
    main()
