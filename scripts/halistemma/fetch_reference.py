#!/usr/bin/env python3
"""Retrieve the pinned host assembly and validate its NCBI sequence inventory."""

import argparse
import gzip
import hashlib
import json
import re
import urllib.request
from pathlib import Path


def sequence_lengths(path: Path) -> dict[str, int]:
    """Validate nonempty, unique nucleotide records without loading chromosomes."""
    lengths: dict[str, int] = {}
    current = None
    with path.open("rb") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                name = line[1:].split()[0].decode()
                if name in lengths:
                    raise ValueError(f"Duplicate sequence ID: {name}")
                current = name
                lengths[name] = 0
            elif current is None or re.search(
                b"[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]", line
            ):
                raise ValueError(f"Invalid nucleotide sequence in {path}: {current}")
            else:
                lengths[current] += len(line)
    if not lengths or not all(lengths.values()):
        raise ValueError(f"Empty reference sequence: {path}")
    return lengths


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--provenance", required=True, type=Path)
    args = parser.parse_args()
    reference = json.loads(args.config.read_text())["host_reference"]
    report = json.loads(Path(reference["dataset_report"]).read_text())["reports"][0]
    if (report["accession"], report["organism"]["tax_id"]) != (
        reference["accession"],
        316183,
    ):
        raise ValueError(
            "Reference accession or organism does not match Halistemma rubrum"
        )
    expected = {}
    mitochondrial = []
    for line in Path(reference["assembly_report"]).read_text().splitlines():
        if line.startswith("#") or not line:
            continue
        fields = line.split("\t")
        expected[fields[4]] = int(fields[8])
        if fields[3] == "Mitochondrion":
            mitochondrial.append(fields[4])
    if not mitochondrial:
        raise ValueError("NCBI sequence inventory has no mitochondrial sequence")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    archive = args.output.with_suffix(".download.gz")
    temporary = args.output.with_suffix(".partial.fna")
    md5 = hashlib.md5()
    with (
        urllib.request.urlopen(reference["url"], timeout=300) as source,
        archive.open("wb") as target,
    ):
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            target.write(block)
            md5.update(block)
    if md5.hexdigest() != reference["md5"]:
        raise ValueError("Downloaded reference does not match NCBI's MD5")
    sha = hashlib.sha256()
    with gzip.open(archive, "rb") as source, temporary.open("wb") as target:
        for block in iter(lambda: source.read(8 * 1024 * 1024), b""):
            target.write(block)
            sha.update(block)
    lengths = sequence_lengths(temporary)
    if lengths != expected:
        raise ValueError(
            "FASTA sequences or lengths disagree with the NCBI assembly report"
        )
    # NCBI reports primary-assembly length separately from the mitochondrial unit.
    primary_bases = sum(lengths.values()) - sum(lengths[name] for name in mitochondrial)
    if primary_bases != int(report["assembly_stats"]["total_sequence_length"]):
        raise ValueError("FASTA length disagrees with the NCBI assembly statistics")
    temporary.replace(args.output)
    args.provenance.write_text(
        json.dumps(
            {
                "status": "PASS",
                "accession": reference["accession"],
                "url": reference["url"],
                "ncbi_md5": md5.hexdigest(),
                "fasta_sha256": sha.hexdigest(),
                "sequences": len(lengths),
                "total_bases": sum(lengths.values()),
                "mitochondrial_sequences": mitochondrial,
            },
            indent=2,
        )
        + "\n"
    )
    archive.unlink()


if __name__ == "__main__":
    main()
