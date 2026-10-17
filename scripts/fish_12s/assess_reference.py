#!/usr/bin/env python3
"""Assess the pilot motif in annotated 12S from long MitoFish assemblies."""

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

MOTIF = "ACAAAC..GGATT..ATACCCC..TATGC"
GROUPS = {
    "117565": "hagfishes",
    "117569": "lampreys",
    "7777": "cartilaginous_fishes",
    "7898": "ray_finned_fishes",
    "118072": "coelacanths",
    "7878": "lungfishes",
}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def fasta_records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        name, chunks = None, []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks).upper()


def closest_window(sequence: str) -> tuple[int, int, str, str]:
    """Count substitutions at fixed motif positions; indels are not allowed."""
    positions = [i for i, base in enumerate(MOTIF) if base != "."]
    expected = np.frombuffer(MOTIF.encode(), dtype=np.uint8)[positions]
    best = (len(positions) + 1, -1, "+", "")
    for strand, seq in [("+", sequence), ("-", reverse_complement(sequence))]:
        if len(seq) < len(MOTIF):
            continue
        bases = np.frombuffer(seq.encode(), dtype=np.uint8)
        windows = np.lib.stride_tricks.sliding_window_view(bases, len(MOTIF))
        distances = np.count_nonzero(windows[:, positions] != expected, axis=1)
        start = int(np.argmin(distances))
        candidate = (
            int(distances[start]),
            start,
            strand,
            seq[start : start + len(MOTIF)],
        )
        if candidate[0] < best[0]:
            best = candidate
    return best


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root, out = args.reference_dir, args.output_dir
    out.mkdir(parents=True, exist_ok=False)
    started = datetime.now(timezone.utc).isoformat()
    annotations = pq.read_table(
        root / "seq_annotation.parquet", filters=[("gene", "=", "12S_rRNA")]
    ).to_pylist()
    # Short amplicons may lack the motif site; do not count them as motif failures.
    annotations = [row for row in annotations if row["length"] >= 700]
    by_accession = defaultdict(list)
    for row in annotations:
        by_accession[row["accession"]].append(row)
    wanted = set(by_accession)
    taxon = {
        r["accession"]: r["taxon_id"]
        for r in pq.read_table(root / "seq_taxonid.parquet").to_pylist()
        if r["accession"] in wanted
    }
    species = {
        r["taxon_id"]: r["species_id"]
        for r in pq.read_table(root / "taxonid_speciesid.parquet").to_pylist()
    }
    names = {
        r["lineage_id"]: r["lineage_name"]
        for r in pq.read_table(root / "taxonid_name.parquet").to_pylist()
    }
    groups, orders = {}, {}
    for row in pq.read_table(root / "speciesid_lineageid.parquet").to_pylist():
        if row["lineage_id"] in GROUPS:
            groups[row["species_id"]] = GROUPS[row["lineage_id"]]
        if row["rank_name"] == "order":
            orders[row["species_id"]] = names.get(row["lineage_id"], row["lineage_id"])
    rows, rejected, seen = [], Counter(), set()
    cache = {}
    with (out / "fish_12s.fasta").open("w") as fasta:
        for header, genome in fasta_records(root / "mitofishdb.fa.gz"):
            if len(genome) < 14000:
                continue
            for accession in header.split()[0].split("|"):
                for annotation in by_accession.get(accession, []):
                    start, end = annotation["start"], annotation["end"]
                    if not (1 <= start <= end <= len(genome)):
                        rejected["invalid_coordinates"] += 1
                        continue
                    seq = genome[start - 1 : end]
                    if annotation["strand"] == "-1":
                        seq = reverse_complement(seq)
                    if len(seq) != annotation["length"]:
                        raise ValueError(f"Coordinate/length disagreement: {accession}")
                    if set(seq) - set("ACGT"):
                        rejected["ambiguous_12s"] += 1
                        continue
                    sp = species.get(taxon.get(accession))
                    if sp not in groups:
                        rejected["unresolved_group"] += 1
                        continue
                    key = (sp, seq)
                    # Equal sequences from one species should not inflate sensitivity.
                    if key in seen:
                        rejected["duplicate_species_sequence"] += 1
                        continue
                    seen.add(key)
                    if seq not in cache:
                        cache[seq] = closest_window(seq)
                    distance, position, strand, window = cache[seq]
                    rid = f"{accession}_{start}_{end}"
                    fasta.write(f">{rid}\n{seq}\n")
                    rows.append(
                        {
                            "reference_id": rid,
                            "accession": accession,
                            "species_id": sp,
                            "species": names[sp],
                            "group": groups[sp],
                            "order": orders.get(sp, "unresolved"),
                            "genome_length": len(genome),
                            "length_12s": len(seq),
                            "min_substitutions": distance,
                            "exact_match": int(distance == 0),
                            "within_one_substitution": int(distance <= 1),
                            "closest_start_0based": position,
                            "closest_strand": strand,
                            "closest_window": window,
                        }
                    )
    write_tsv(out / "reference_matches.tsv", rows)
    summary = []
    for level in ["group", "order"]:
        values = sorted({r[level] for r in rows})
        if level == "group":
            values = ["all"] + values
        for value in values:
            subset = rows if value == "all" else [r for r in rows if r[level] == value]
            summary.append(
                {
                    "level": level,
                    "name": value,
                    "reference_sequences": len(subset),
                    "species": len({r["species_id"] for r in subset}),
                    "exact_sequences": sum(r["exact_match"] for r in subset),
                    "one_substitution_sequences": sum(
                        r["within_one_substitution"] for r in subset
                    ),
                    "species_with_exact": len(
                        {r["species_id"] for r in subset if r["exact_match"]}
                    ),
                    "species_with_one_substitution": len(
                        {
                            r["species_id"]
                            for r in subset
                            if r["within_one_substitution"]
                        }
                    ),
                }
            )
    write_tsv(out / "coverage_summary.tsv", summary)
    inputs = [root / "mitofishdb.fa.gz", *sorted(root.glob("*.parquet"))]
    inputs = [
        p
        for p in inputs
        if p.name not in {"seq_heterospecific.parquet", "seq_description.parquet"}
    ]
    record = {
        "started_utc": started,
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "command": sys.argv,
        "python": sys.version,
        "numpy": np.__version__,
        "pyarrow": pq.__version__
        if hasattr(pq, "__version__")
        else __import__("pyarrow").__version__,
        "script_sha256": sha256(Path(__file__)),
        "motif": MOTIF,
        "reference_release": "MitoFish 2026.09",
        "selection": "genome >=14000 nt; annotated 12S >=700 nt; ACGT only; unique sequence per species",
        "excluded": dict(rejected),
        "inputs": {str(p): sha256(p) for p in inputs},
        "outputs": {str(p): sha256(p) for p in out.iterdir() if p.is_file()},
    }
    (out / "run.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(summary[:7], indent=2))


if __name__ == "__main__":
    main()
