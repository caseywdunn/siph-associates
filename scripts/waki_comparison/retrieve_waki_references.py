#!/usr/bin/env python3
"""Archive Waki et al. 2026 marker sequences and their source annotations.

Requires Python >=3.10, Biopython and curl. Existing downloads are reused;
delete the source file explicitly to request a fresh copy. No alignments run.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import re
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from zipfile import ZipFile

import Bio
from Bio import SeqIO

TITLE = (
    "Metacercariae infecting seven cnidarian species with their life cycle "
    "information including their adult stages in Japan"
)
DOI = "10.1017/S0022149X25100989"
PUBMED_LINK_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?"
    "dbfrom=pubmed&db=nuccore&id=41424079&retmode=json"
)
SUPPLEMENT_URL = (
    "https://static.cambridge.org/content/id/"
    "urn%3Acambridge.org%3Aid%3Aarticle%3AS0022149X25100989/"
    "resource/name/S0022149X25100989sup001.xlsx"
)
SUPPLEMENT_SHA256 = "e591d82bbbd510d5b7e9df3b8f4695df8965eafafa92b4fcb5ab81fdab902ab8"
# These eight accessions are explicit in the article's DNA-marker lists but
# absent from its supplementary distance matrices. Do not fill numeric gaps.
PAPER_ONLY_ACCESSIONS = {
    "LC889172": ("28S", "Tetrochetus coryphaenae", "p10"),
    "LC889180": ("28S", "Cephalolepidapedon saba", "p15"),
    "LC889181": ("28S", "Cephalolepidapedon saba", "p15"),
    "LC889188": ("28S", "Lepocreadiidae sp. 1", "p18"),
    "LC889268": ("COI", "Opechona sp. 1", "p15"),
    "LC889269": ("COI", "Opechona sp. 1", "p15"),
    "LC889270": ("COI", "Lepocreadiidae sp. 1", "p18"),
    "LC889271": ("COI", "Lepocreadiidae sp. 1", "p18"),
}
# Independently found through the publication's NCBI sequence links. Both
# records cite this article, although its printed accession lists omit them.
DATABASE_ONLY_ACCESSIONS = {
    "LC889189": ("28S", "Lepocreadiidae sp. 1"),
    "LC889237": ("COI", "Dinurus barbatus"),
}
NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def sha256(path: Path) -> str:
    """Return a file digest."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download(url: str, path: Path) -> dict[str, str]:
    """Fetch a missing source atomically or record reuse of its cached bytes."""
    reused = path.exists()
    if not reused:
        temporary = path.with_suffix(path.suffix + ".download")
        subprocess.run(
            [
                "curl",
                "--fail",
                "--location",
                "--max-time",
                "90",
                url,
                "-o",
                str(temporary),
            ],
            check=True,
        )
        temporary.replace(path)
    return {
        "url": url,
        "path": str(path),
        "sha256": sha256(path),
        "cache_reused": str(reused).lower(),
        "retrieved_or_cached_file_mtime_utc": datetime.fromtimestamp(
            path.stat().st_mtime, timezone.utc
        ).isoformat(),
    }


def supplement_rows(path: Path) -> tuple[list[dict], list[dict]]:
    """Read accession rows from the archived XLSX using its cached cell values."""
    mapped, unmapped = [], []
    with ZipFile(path) as archive:
        strings = [
            "".join(item.itertext())
            for item in ET.fromstring(archive.read("xl/sharedStrings.xml"))
        ]
        links = {
            element.attrib["Id"]: element.attrib["Target"]
            for element in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        }
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        for sheet in workbook.find("s:sheets", NS):
            target = links[sheet.attrib[f"{{{REL_NS}}}id"]]
            target = target.lstrip("/") if target.startswith("/") else f"xl/{target}"
            root = ET.fromstring(archive.read(target))
            marker = ""
            for row in root.findall("s:sheetData/s:row", NS):
                cells = {}
                for cell in row.findall("s:c", NS):
                    value = cell.find("s:v", NS)
                    if value is None or value.text is None:
                        continue
                    text = (
                        strings[int(value.text)] if cell.get("t") == "s" else value.text
                    )
                    cells[re.sub(r"\d", "", cell.attrib["r"])] = text.strip()
                if "Supplementary Table" in cells.get("A", ""):
                    title = cells["A"]
                    marker = next(x for x in ("28S", "ITS2", "COI") if x in title)
                label = cells.get("B", "")
                if not label or not cells.get("A", "").isdigit():
                    continue
                if not marker:
                    raise ValueError(f"Missing marker for {sheet.attrib['name']}")
                accession_text = cells.get("C", "")
                base = {
                    "sheet": sheet.attrib["name"].strip(),
                    "xlsx_row": row.attrib["r"],
                    "sequence_label": label,
                    "marker": marker,
                    "original_accession_cell": accession_text,
                }
                accessions = re.findall(
                    r"\b[A-Z]{1,3}\d{5,9}(?:\.\d+)?\b", accession_text
                )
                if not accessions:
                    unmapped.append(base)
                for accession in accessions:
                    mapped.append({**base, "accession": accession})
    return mapped, unmapped


def write_tsv(path: Path, rows: list[dict]) -> None:
    """Write a nonempty, ordered rectangular table."""
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("data/sources/waki2026")
    )
    args = parser.parse_args()
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    supplementary = output / "S0022149X25100989sup001.xlsx"
    sources = [download(SUPPLEMENT_URL, supplementary)]
    if sha256(supplementary) != SUPPLEMENT_SHA256:
        raise ValueError(
            "Supplement changed; inspect its accession assignments before proceeding"
        )
    mapped, unmapped = supplement_rows(supplementary)
    write_tsv(output / "supplement_accessions.tsv", mapped)
    write_tsv(output / "supplement_unaccessioned.tsv", unmapped)
    rows_by_accession = defaultdict(list)
    for row in mapped:
        if row["accession"].startswith("LC889"):
            rows_by_accession[row["accession"].split(".")[0]].append(row)
    accessions = sorted(
        set(rows_by_accession)
        | set(PAPER_ONLY_ACCESSIONS)
        | set(DATABASE_ONLY_ACCESSIONS)
    )
    if len(accessions) != 108:
        raise ValueError(
            f"Expected 108 verified study accessions; found {len(accessions)}"
        )
    links_path = output / "pubmed_nuccore_links.json"
    sources.append(download(PUBMED_LINK_URL, links_path))
    links = json.loads(links_path.read_text())
    sequence_ids = [
        uid
        for linkset in links["linksets"]
        for database in linkset.get("linksetdbs", [])
        if database["linkname"] == "pubmed_nuccore"
        for uid in database["links"]
    ]
    if len(sequence_ids) != 108 or len(set(sequence_ids)) != 108:
        raise ValueError(
            "The publication's NCBI sequence links changed; review accession completeness"
        )
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(
        {
            "db": "nuccore",
            "id": ",".join(sequence_ids),
            "rettype": "gb",
            "retmode": "text",
        }
    )
    genbank = output / "references.gb"
    sources.append(download(url, genbank))
    records = sorted(SeqIO.parse(genbank, "genbank"), key=lambda record: record.id)
    if len(records) != len(accessions) or {r.id.split(".")[0] for r in records} != set(
        accessions
    ):
        raise ValueError("Retrieved accessions do not match requested accession set")
    if len({r.id for r in records}) != len(records):
        raise ValueError("Duplicate versioned sequence identifiers")
    metadata, annotations = [], []
    for record in records:
        accession = record.id.split(".")[0]
        titles = [ref.title for ref in record.annotations.get("references", [])]
        if TITLE not in titles:
            raise ValueError(f"{record.id} does not cite the Waki article")
        sequence = str(record.seq).upper()
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Invalid nucleotide sequence: {record.id}")
        source_features = [
            feature for feature in record.features if feature.type == "source"
        ]
        if len(source_features) != 1:
            raise ValueError(f"Unexpected source-feature count: {record.id}")
        qualifiers = source_features[0].qualifiers
        marker = (
            "ITS2"
            if "ITS2" in record.description
            else "28S"
            if "28S" in record.description
            else "COI"
        )
        assignments = rows_by_accession.get(accession, [])
        if any(row["marker"] != marker for row in assignments):
            raise ValueError(f"Marker mismatch for {record.id}")
        paper_only = PAPER_ONLY_ACCESSIONS.get(accession)
        database_only = DATABASE_ONLY_ACCESSIONS.get(accession)
        if paper_only and marker != paper_only[0]:
            raise ValueError(f"Article marker mismatch for {record.id}")
        if database_only and marker != database_only[0]:
            raise ValueError(f"Database-only marker mismatch for {record.id}")
        labels = sorted({row["sequence_label"] for row in assignments})
        taxon = (
            re.split(r"\s*\(", labels[0], maxsplit=1)[0].strip()
            if labels
            else (paper_only or database_only)[1]
        )
        metadata.append(
            {
                "accession_version": record.id,
                "accession": accession,
                "marker": marker,
                "paper_taxon": taxon,
                "genbank_organism": record.annotations["organism"],
                "host": "|".join(qualifiers.get("host", [])),
                "isolate": "|".join(qualifiers.get("isolate", [])),
                "developmental_stage": "|".join(qualifiers.get("dev_stage", [])),
                "location": "|".join(qualifiers.get("geo_loc_name", [])),
                "collection_date": "|".join(qualifiers.get("collection_date", [])),
                "length_bp": len(record),
                "supplement_sheets": "|".join(
                    sorted({row["sheet"] for row in assignments})
                ),
                "supplement_labels": "|".join(labels),
                "paper_only_accession_location": paper_only[2] if paper_only else "",
                "database_only_accession_evidence": "NCBI pubmed_nuccore link for PMID41424079; GenBank article title verified"
                if database_only
                else "",
                "description": record.description,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
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
    write_tsv(output / "reference_metadata.tsv", metadata)
    SeqIO.write(records, output / "references.fasta", "fasta")
    (output / "reference_features.json").write_text(
        json.dumps(annotations, indent=2) + "\n"
    )
    artifacts = [
        supplementary,
        links_path,
        genbank,
        output / "references.fasta",
        output / "reference_metadata.tsv",
        output / "reference_features.json",
        output / "supplement_accessions.tsv",
        output / "supplement_unaccessioned.tsv",
    ]
    manifest = {
        "status": "complete",
        "publication_doi": DOI,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "command": shlex.join([sys.executable, *sys.argv]),
        "script_sha256": sha256(Path(__file__)),
        "python": platform.python_version(),
        "biopython": Bio.__version__,
        "sources": sources,
        "selection": "108 NCBI sequences linked to PMID41424079:98 supplement accessions,eight additional printed article accessions,and two database-only accessions; every record cites the article title",
        "record_count": len(records),
        "marker_counts": dict(Counter(row["marker"] for row in metadata)),
        "host_counts": dict(Counter(row["host"] for row in metadata)),
        "output_sha256": {path.name: sha256(path) for path in artifacts},
        "limitations": [
            "No 18S sequences are deposited in this study reference set.",
            "ITS2 records span 5.8S/ITS2/28S; their internal boundaries are not annotated.",
            "Unaccessioned supplement rows are preserved separately, not invented as sequences.",
            "Specimen-specific GenBank hosts and life stages are preserved without inferring other hosts.",
        ],
    }
    (output / "retrieval_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(json.dumps({"records": len(records), "markers": manifest["marker_counts"]}))


if __name__ == "__main__":
    main()
