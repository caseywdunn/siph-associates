#!/usr/bin/env python3
"""Prepare validated MAG sequences, GTDB genes, and broad reference metadata."""

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
from pathlib import Path

from Bio import SeqIO


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_parents(path: Path) -> dict[int, int]:
    parents = {}
    with path.open() as handle:
        for line in handle:
            fields = line.split("|")
            parents[int(fields[0])] = int(fields[1])
    return parents


def domain(taxid: int, parents: dict[int, int]) -> str:
    domains = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota", 10239: "Viruses"}
    visited = set()
    while taxid not in visited:
        if taxid in domains:
            return domains[taxid]
        visited.add(taxid)
        taxid = parents.get(taxid, taxid)
    return "unresolved_taxonomy"


def prepare_inputs(config: dict, root: Path) -> None:
    destination = root / "inputs"
    destination.mkdir(parents=True, exist_ok=True)
    contigs, tiles = [], []
    seen = set()
    with (root / "tiles.fna").open("w") as output:
        for item in config["bins"]:
            source = Path(item["fasta"])
            if hashlib.sha256(source.read_bytes()).hexdigest() != item["sha256"]:
                raise ValueError(f"Input checksum changed: {source}")
            records = list(SeqIO.parse(source, "fasta"))
            if not records:
                raise ValueError(f"Empty FASTA: {source}")
            for record in records:
                sequence = str(record.seq).upper()
                if record.id in seen or not sequence or set(sequence) - set("ACGTN"):
                    raise ValueError(f"Invalid or duplicate contig: {record.id}")
                seen.add(record.id)
                contigs.append(
                    {
                        "contig_id": record.id,
                        "bin_id": item["bin_id"],
                        "length_nt": len(sequence),
                    }
                )
                for start in range(0, len(sequence), config["tile_step"]):
                    end = min(start + config["tile_length"], len(sequence))
                    if end - start < 150:
                        break
                    query = f"{item['bin_id']}@{record.id}:{start + 1}-{end}"
                    output.write(f">{query}\n{sequence[start:end]}\n")
                    tiles.append(
                        {
                            "query_id": query,
                            "bin_id": item["bin_id"],
                            "contig_id": record.id,
                            "start": start + 1,
                            "end": end,
                            "strand": "+",
                            "length": end - start,
                        }
                    )
            shutil.copyfile(source, destination / f"{item['bin_id']}.fa")
    write_tsv(root / "contigs.tsv", contigs, ["contig_id", "bin_id", "length_nt"])
    write_tsv(
        root / "tiles.tsv",
        tiles,
        ["query_id", "bin_id", "contig_id", "start", "end", "strand", "length"],
    )
    (root / "inputs.json").write_text(json.dumps(config["bins"], indent=2) + "\n")


def prepare_genes(config: dict, root: Path) -> None:
    genes, markers = [], []
    seen = set()
    marker_root = root / "gtdb_identify/identify/intermediate_results/marker_genes"
    with (root / "genes.faa").open("w") as output:
        for item in config["bins"]:
            bin_id = item["bin_id"]
            folder = marker_root / bin_id
            ids = {}
            for record in SeqIO.parse(folder / f"{bin_id}_protein.faa", "fasta"):
                query = f"{bin_id}@{record.id}"
                if query in seen:
                    raise ValueError(f"Duplicate gene: {query}")
                seen.add(query)
                ids[record.id] = query
                fields = record.description.split(" # ")
                sequence = str(record.seq).rstrip("*")
                if not sequence or "*" in sequence:
                    raise ValueError(f"Invalid predicted protein: {query}")
                output.write(f">{query}\n{sequence}\n")
                genes.append(
                    {
                        "query_id": query,
                        "bin_id": bin_id,
                        "contig_id": record.id.rsplit("_", 1)[0],
                        "start": int(fields[1]),
                        "end": int(fields[2]),
                        "strand": fields[3],
                        "length": len(sequence),
                    }
                )
            if not ids:
                raise ValueError(f"No genes produced: {bin_id}")
            for suffix in ["pfam", "tigrfam"]:
                with (folder / f"{bin_id}_{suffix}_tophit.tsv").open() as handle:
                    next(handle)
                    for line in handle:
                        gene, hits = line.rstrip().split("\t")
                        for hit in hits.split(";"):
                            marker, evalue, bitscore = hit.split(",")
                            if marker in config["archaeal_marker_ids"]:
                                markers.append(
                                    {
                                        "query_id": ids[gene],
                                        "bin_id": bin_id,
                                        "marker_id": marker,
                                        "hmm_evalue": evalue,
                                        "hmm_bitscore": bitscore,
                                    }
                                )
    write_tsv(
        root / "genes.tsv",
        genes,
        ["query_id", "bin_id", "contig_id", "start", "end", "strand", "length"],
    )
    write_tsv(
        root / "archaeal_marker_hits.tsv",
        markers,
        ["query_id", "bin_id", "marker_id", "hmm_evalue", "hmm_bitscore"],
    )


def prepare_reference(root: Path) -> None:
    ref = root / "reference"
    parents = read_parents(ref / "taxonomy/nodes.dmp")
    rows, controls = [], []
    counts = {"Bacteria": 0, "Archaea": 0, "Eukaryota": 0}
    with (
        gzip.open(ref / "uniprot_sprot.fasta.gz", "rt") as handle,
        (ref / "controls.faa").open("w") as output,
    ):
        for record in SeqIO.parse(handle, "fasta"):
            match = re.search(r" OX=(\d+)", record.description)
            if not match:
                raise ValueError(f"Swiss-Prot entry lacks taxid: {record.id}")
            taxid = int(match[1])
            organism = re.search(r" OS=(.+?) OX=", record.description).group(1)
            label = domain(taxid, parents)
            rows.append(
                {
                    "subject_id": record.id,
                    "accession": record.id.split("|")[1],
                    "taxid": taxid,
                    "domain": label,
                    "organism": organism,
                    "description": record.description,
                }
            )
            if label in counts and counts[label] < 5 and 100 <= len(record) <= 500:
                query = f"CONTROL_{label}_{counts[label]}"
                output.write(f">{query}\n{record.seq}\n")
                controls.append(
                    {
                        "query_id": query,
                        "bin_id": "control",
                        "contig_id": record.id,
                        "start": 1,
                        "end": len(record),
                        "strand": "+",
                        "length": len(record),
                        "expected_domain": label,
                    }
                )
                counts[label] += 1
    if set(counts.values()) != {5}:
        raise ValueError(f"Missing reference controls: {counts}")
    write_tsv(
        ref / "swissprot_metadata.tsv",
        rows,
        ["subject_id", "accession", "taxid", "domain", "organism", "description"],
    )
    write_tsv(
        ref / "controls.tsv",
        controls,
        [
            "query_id",
            "bin_id",
            "contig_id",
            "start",
            "end",
            "strand",
            "length",
            "expected_domain",
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["inputs", "genes", "reference"])
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    root = Path(config["run"])
    root.mkdir(parents=True, exist_ok=True)
    if args.stage == "inputs":
        prepare_inputs(config, root)
    elif args.stage == "genes":
        prepare_genes(config, root)
    else:
        prepare_reference(root)


if __name__ == "__main__":
    main()
