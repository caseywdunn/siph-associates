#!/usr/bin/env python3
"""Archive the external comparison accessions in Waki et al. 2026 supplement.

Run retrieve_waki_references.py first. Keeps external references separate from
the 106 sequences newly deposited for the Waki study; performs no alignment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import shlex
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

import Bio
from Bio import SeqIO
from retrieve_waki_references import download, sha256, write_tsv


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/sources/waki2026")
    )
    args = parser.parse_args()
    output = args.output_dir
    assignments = output / "supplement_accessions.tsv"
    rows_by_accession = defaultdict(list)
    with assignments.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if not row["accession"].startswith("LC889"):
                rows_by_accession[row["accession"].split(".")[0]].append(row)
    accessions = sorted(rows_by_accession)
    if len(accessions) != 113:
        raise ValueError(f"Expected 113 external accessions; found {len(accessions)}")
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(
        {
            "db": "nuccore",
            "id": ",".join(accessions),
            "rettype": "gb",
            "retmode": "text",
        }
    )
    genbank = output / "context_references.gb"
    source = download(url, genbank)
    records = list(SeqIO.parse(genbank, "genbank"))
    if len(records) != len(accessions) or {r.id.split(".")[0] for r in records} != set(
        accessions
    ):
        raise ValueError("Retrieved context accessions do not match requested set")
    if len({r.id for r in records}) != len(records):
        raise ValueError("Duplicate versioned context identifiers")
    metadata, annotations = [], []
    for record in records:
        accession = record.id.split(".")[0]
        rows = rows_by_accession[accession]
        markers = {row["marker"] for row in rows}
        if len(markers) != 1:
            raise ValueError(
                f"Conflicting supplement markers for {record.id}: {markers}"
            )
        marker = next(iter(markers))
        sequence = str(record.seq).upper()
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Invalid nucleotide sequence: {record.id}")
        sources = [feature for feature in record.features if feature.type == "source"]
        if len(sources) != 1:
            raise ValueError(f"Unexpected source-feature count for {record.id}")
        qualifiers = sources[0].qualifiers
        labels = sorted({row["sequence_label"] for row in rows})
        paper_taxon = re.split(r"\s*\(", labels[0], maxsplit=1)[0].strip()
        metadata.append(
            {
                "accession_version": record.id,
                "accession": accession,
                "marker": marker,
                "paper_taxon": paper_taxon,
                "genbank_organism": record.annotations["organism"],
                "host": "|".join(qualifiers.get("host", [])),
                "isolate": "|".join(qualifiers.get("isolate", [])),
                "developmental_stage": "|".join(qualifiers.get("dev_stage", [])),
                "location": "|".join(qualifiers.get("geo_loc_name", [])),
                "collection_date": "|".join(qualifiers.get("collection_date", [])),
                "length_bp": len(record),
                "supplement_sheets": "|".join(sorted({row["sheet"] for row in rows})),
                "supplement_labels": "|".join(labels),
                "description": record.description,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "reference_titles": "|".join(
                    ref.title for ref in record.annotations.get("references", [])
                ),
            }
        )
        for feature in record.features:
            annotations.append(
                {
                    "accession_version": record.id,
                    "type": feature.type,
                    "location_biopython": str(feature.location),
                    "start_zero_based": int(feature.location.start),
                    "end_exclusive": int(feature.location.end),
                    "strand": feature.location.strand,
                    "qualifiers": feature.qualifiers,
                }
            )
    write_tsv(output / "context_metadata.tsv", metadata)
    SeqIO.write(records, output / "context_references.fasta", "fasta")
    (output / "context_features.json").write_text(
        json.dumps(annotations, indent=2) + "\n"
    )
    artifacts = [
        genbank,
        output / "context_references.fasta",
        output / "context_metadata.tsv",
        output / "context_features.json",
    ]
    manifest = {
        "status": "complete",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command": shlex.join([sys.executable, *sys.argv]),
        "scripts_sha256": {
            Path(__file__).name: sha256(Path(__file__)),
            "retrieve_waki_references.py": sha256(
                Path(__file__).with_name("retrieve_waki_references.py")
            ),
        },
        "input_sha256": {assignments.name: sha256(assignments)},
        "source": source,
        "python": platform.python_version(),
        "biopython": Bio.__version__,
        "record_count": len(records),
        "marker_counts": dict(Counter(row["marker"] for row in metadata)),
        "output_sha256": {path.name: sha256(path) for path in artifacts},
        "limitations": [
            "External supplement references are distinct from newly deposited Waki study specimens.",
            "Supplement marker labels may name only part of a longer deposited amplicon.",
            "GenBank host, isolate or other qualifiers may be absent; missing values are not inferred.",
            "Paper labels are preserved alongside current GenBank names without resolving taxonomic synonyms.",
        ],
    }
    (output / "context_retrieval_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps({"records": len(records), "markers": manifest["marker_counts"]}))


if __name__ == "__main__":
    main()
