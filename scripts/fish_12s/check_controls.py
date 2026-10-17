#!/usr/bin/env python3
"""Verify reference extraction, non-fish specificity, and paired filter behavior."""

import argparse
import csv
import json
import random
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

from pilot import MOTIF, digest, pair_ids, read_id, records


def write_table(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--assessment-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bbduk", required=True)
    parser.add_argument("--search-reference", type=Path)
    args = parser.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=False)
    reference = args.assessment_dir / "fish_12s.fasta"
    refs = {r.id: str(r.seq) for r in SeqIO.parse(reference, "fasta")}
    rows = []
    for rec in SeqIO.parse(args.reference_dir / "fish_spotchecks.gb", "genbank"):
        for rid, seq in refs.items():
            if not rid.startswith(rec.id.split(".")[0] + "_"):
                continue
            genome = str(rec.seq) * 2  # MitoFish may rotate circular genomes.
            valid = seq in genome or str(Seq(seq).reverse_complement()) in genome
            rows.append(
                {
                    "accession": rec.id,
                    "organism": rec.annotations["organism"],
                    "reference_id": rid,
                    "matches_original_genome": valid,
                }
            )
    if len(rows) != 6 or not all(r["matches_original_genome"] for r in rows):
        raise ValueError("GenBank cross-check failed")
    write_table(out / "genbank_spotchecks.tsv", rows)
    rows, frog = [], None
    for rec in SeqIO.parse(args.reference_dir / "nonfish_controls.gb", "genbank"):
        hits = list(re.finditer(MOTIF, str(rec.seq)))
        rows.append(
            {
                "accession": rec.id,
                "organism": rec.annotations["organism"],
                "genome_length": len(rec),
                "exact_motif_matches": len(hits),
                "matched_sequences": ";".join(m.group() for m in hits),
            }
        )
        if rec.id == "NC_001573.1":
            if len(hits) != 1:
                raise ValueError("Frog positive control changed")
            frog = str(rec.seq)[hits[0].start() - 50 : hits[0].start() + 101]
    write_table(out / "nonfish_controls.tsv", rows)
    with (args.assessment_dir / "reference_matches.tsv").open() as handle:
        assessed = list(csv.DictReader(handle, delimiter="\t"))
    ray = next(
        r
        for r in assessed
        if r["group"] == "ray_finned_fishes" and r["exact_match"] == "1"
    )
    shark = next(r for r in assessed if r["group"] == "cartilaginous_fishes")
    rayseq, sharkseq = refs[ray["reference_id"]], refs[shark["reference_id"]]
    start = re.search(MOTIF, rayseq).start()
    rayread = rayseq[start - 50 : start + 101]
    sharkstart = int(shark["closest_start_0based"])
    sharkread = sharkseq[sharkstart - 50 : sharkstart + 101]
    rng = random.Random(20261017)
    random_read = "".join(rng.choices("ACGT", k=151))
    quality_motif = re.search(MOTIF, rayread).group()
    quality = "I" * 50 + quality_motif + "I" * (101 - len(quality_motif))
    examples = [
        ("motif_R1", rayread, random_read, "I" * 151, "I" * 151),
        (
            "motif_reverse_R2",
            random_read,
            str(Seq(rayread).reverse_complement()),
            "I" * 151,
            "I" * 151,
        ),
        ("fish_away_from_motif", rayseq[:151], random_read, "I" * 151, "I" * 151),
        ("cartilaginous_motif_variant", sharkread, random_read, "I" * 151, "I" * 151),
        ("frog", frog, random_read, "I" * 151, "I" * 151),
        ("quality_only", random_read, random_read, quality, "I" * 151),
        ("negative", random_read, random_read, "I" * 151, "I" * 151),
    ]
    if args.search_reference:
        examples.extend(
            [
                ("poly_G", "G" * 151, "G" * 151, "I" * 151, "I" * 151),
                (
                    "ATTT_repeat",
                    ("ATTT" * 38)[:151],
                    ("AAAT" * 38)[:151],
                    "I" * 151,
                    "I" * 151,
                ),
            ]
        )
    for mate in [1, 2]:
        with (out / f"input_R{mate}.fastq").open("w") as handle:
            for example in examples:
                name, seq, qual = example[0], example[mate], example[mate + 2]
                if len(seq) != 151 or len(qual) != 151:
                    raise ValueError("Invalid synthetic read length")
                handle.write(f"@{name}/{mate}\n{seq}\n+\n{qual}\n")
    commands = []
    for mate in [1, 2]:
        cmd = [
            "seqkit",
            "grep",
            "-j",
            "1",
            "-s",
            "-r",
            "-P",
            "-p",
            MOTIF,
            str(out / f"input_R{mate}.fastq"),
        ]
        with (out / f"motif_R{mate}.fastq").open("w") as handle:
            subprocess.run(cmd, stdout=handle, check=True)
        commands.append(cmd)
    cmd = [
        args.bbduk,
        "-Xmx1g",
        "threads=1",
        "k=31",
        "hdist=0",
        "maskmiddle=f",
        "rcomp=t",
        "rieb=t",
        "qtrim=f",
        "minlen=0",
        "ordered=t",
        "overwrite=f",
        f"ref={args.search_reference or reference}",
        f"in={out / 'input_R1.fastq'}",
        f"in2={out / 'input_R2.fastq'}",
        f"outm={out / 'kmer_R1.fastq'}",
        f"outm2={out / 'kmer_R2.fastq'}",
    ]
    with (out / "bbduk.log").open("w") as handle:
        subprocess.run(cmd, stdout=handle, stderr=handle, check=True)
    commands.append(cmd)
    motif_ids = {
        read_id(name)
        for mate in [1, 2]
        for name, _, _ in records(out / f"motif_R{mate}.fastq")
    }
    if motif_ids != {"motif_R1", "motif_reverse_R2", "frog"}:
        raise ValueError(f"Unexpected regex results: {motif_ids}")
    kmer_ids = pair_ids(out / "kmer_R1.fastq", out / "kmer_R2.fastq")
    fish_ids = {
        "motif_R1",
        "motif_reverse_R2",
        "fish_away_from_motif",
        "cartilaginous_motif_variant",
    }
    if not fish_ids <= kmer_ids or kmer_ids & {"negative", "quality_only"}:
        raise ValueError(f"Unexpected kmer results: {kmer_ids}")
    if args.search_reference and kmer_ids & {"poly_G", "ATTT_repeat"}:
        raise ValueError("Masked reference retained a repetitive control")
    rows = [
        {
            "control": r[0],
            "regex_pair": int(r[0] in motif_ids),
            "whole_12s_pair": int(r[0] in kmer_ids),
        }
        for r in examples
    ]
    write_table(out / "filter_controls.tsv", rows)
    provenance = {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "commands": commands,
        "script_sha256": digest(Path(__file__)),
        "helper_sha256": digest(Path(__file__).with_name("pilot.py")),
        "inputs": {
            str(p): digest(p)
            for p in [
                reference,
                args.assessment_dir / "reference_matches.tsv",
                args.reference_dir / "nonfish_controls.gb",
                args.reference_dir / "fish_spotchecks.gb",
            ]
        },
        "outputs": {str(p): digest(p) for p in out.iterdir() if p.is_file()},
    }
    if args.search_reference:
        provenance["inputs"][str(args.search_reference)] = digest(args.search_reference)
    (out / "run.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
