#!/usr/bin/env python3
"""Annotate the frozen competitive references from accession-linked records."""

import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

from Bio import SeqIO

from pilot import digest

FISH_LINEAGES = {
    "Actinopterygii",
    "Actinopteri",
    "Chondrichthyes",
    "Myxini",
    "Petromyzontida",
    "Petromyzonti",
    "Petromyzontiformes",
    "Dipnoi",
    "Dipnomorpha",
    "Coelacanthimorpha",
    "Actinistia",
}


def is_fish(lineage: list[str]) -> bool:
    return bool(FISH_LINEAGES.intersection(lineage))


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--genbank", type=Path, required=True)
    parser.add_argument("--mitofish", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    frozen = {}
    with gzip.open(args.fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            accession = record.id.split("|")[1] if "|" in record.id else record.id
            if accession in frozen:
                raise ValueError(f"Duplicate frozen accession: {accession}")
            frozen[accession] = (record.id, sequence_hash(str(record.seq)), len(record))
    rows, seen, contig_only = [], set(), []
    with gzip.open(args.genbank, "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            if record.id not in frozen:
                continue
            rid, expected, length = frozen[record.id]
            if record.id in seen or len(record) != length:
                raise ValueError(f"GenBank does not match frozen sequence: {record.id}")
            if record.seq.defined:
                if sequence_hash(str(record.seq)) != expected:
                    raise ValueError(f"GenBank sequence changed: {record.id}")
            elif record.annotations.get("contig"):
                # Some RefSeq GenBank records provide CONTIG rather than ORIGIN.
                # Their immutable accession.version and length still identify the record.
                contig_only.append(record.id)
            else:
                raise ValueError(
                    f"GenBank record has neither sequence nor CONTIG: {record.id}"
                )
            seen.add(record.id)
            lineage = record.annotations["taxonomy"]
            taxids = [
                v.removeprefix("taxon:")
                for f in record.features
                if f.type == "source"
                for v in f.qualifiers.get("db_xref", [])
                if v.startswith("taxon:")
            ]
            intervals = [
                f"{int(part.start) + 1}-{int(part.end)}"
                for f in record.features
                if f.type == "rRNA"
                and any(
                    "12S" in value.upper() for value in f.qualifiers.get("product", [])
                )
                for part in f.location.parts
            ]
            rows.append(
                dict(
                    reference_id=rid,
                    accession=record.id,
                    name=record.annotations["organism"],
                    taxid=";".join(taxids),
                    fish=int(is_fish(lineage)),
                    lineage=";".join(lineage),
                    marker_intervals_1based=";".join(intervals),
                    source="RefSeq_GenBank",
                )
            )
    if seen != frozen.keys():
        raise ValueError(
            f"Missing frozen accessions in GenBank: {sorted(frozen.keys() - seen)[:20]}"
        )
    with args.mitofish.open() as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                dict(
                    reference_id=r["reference_id"],
                    accession=r["accession"],
                    name=r["species"],
                    taxid=r["species_id"],
                    fish=1,
                    lineage=r["group"],
                    marker_intervals_1based=f"1-{r['length_12s']}",
                    source="MitoFish_annotated_12S",
                )
            )
    if len({r["reference_id"] for r in rows}) != len(rows):
        raise ValueError("Competitive reference identifiers are not unique")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        writer = csv.DictWriter(
            handle, list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    args.output.with_suffix(".json").write_text(
        json.dumps(
            dict(
                refseq_records=len(seen),
                total_records=len(rows),
                fish_records=sum(r["fish"] for r in rows),
                input_sha256={
                    str(p): digest(p)
                    for p in [args.fasta, args.genbank, args.mitofish, Path(__file__)]
                },
                output_sha256=digest(args.output),
                sequence_hash_checked=len(seen) - len(contig_only),
                accession_version_and_length_checked=contig_only,
            ),
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
