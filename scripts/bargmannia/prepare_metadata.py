#!/usr/bin/env python3
"""Link only the eight new Table S2 libraries to public SRA and local reads."""

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from datetime import UTC, datetime, timedelta
from pathlib import Path

STUDY = "Manko2026"
DOI = "10.1093/zoolinnean/zlag150"
RAW = Path("/gpfs/ycga/work/dunn/sequences/illumina/apolemia_bargmannia_maciej")
SOURCE = Path("data/sources/manko2026")
OUT = Path("data/metadata/manko2026")
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def sheet(data):
    with zipfile.ZipFile(io.BytesIO(data)) as book:
        strings = [
            "".join(e.itertext())
            for e in ET.fromstring(book.read("xl/sharedStrings.xml")).findall(
                "m:si", NS
            )
        ]
        rows = []
        for element in ET.fromstring(book.read("xl/worksheets/sheet1.xml")).findall(
            ".//m:row", NS
        ):
            cells = {}
            for cell in element:
                node = cell.find("m:v", NS)
                value = node.text if node is not None else ""
                if cell.get("t") == "s":
                    value = strings[int(value)]
                cells[re.sub(r"\d+", "", cell.get("r"))] = value.strip()
            rows.append((int(element.get("r")), cells))
    header = rows[0][1]
    return [
        dict(
            source_row=n, **{key: r.get(col, "") for col, key in header.items() if key}
        )
        for n, r in rows[1:]
    ]


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--download", action="store_true", help="Refresh public SRA XML metadata"
    )
    args = parser.parse_args()
    SOURCE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    supplement = Path("tmp/zlag150_supplementary_data.zip")
    with zipfile.ZipFile(supplement) as archive:
        s1 = sheet(archive.read("Table S1.xlsx"))
        s2 = sheet(archive.read("Table S2.xlsx"))
    selected = [
        r
        for r in s2
        if r["Source"] == "This study" and r["Raw reads"].startswith("SRR")
    ]
    assert len(selected) == 8
    assert {r["Raw reads"] for r in selected} == {
        f"SRR{i}" for i in range(36270048, 36270056)
    }
    assert {r["Species"] for r in selected} == {
        "Bargmannia " + x
        for x in ["amoena", "elongata", "marlena", "profunda", "pyrosa", "stenotes"]
    }
    write(OUT / "table_s1.tsv", s1)
    write(OUT / "table_s2.tsv", s2)
    metadata, files, inputs = (
        [],
        [],
        [supplement, Path("tmp/zlag150.pdf"), Path(__file__)],
    )
    for row in selected:
        event, run, species = row["Specimen_ID"], row["Raw reads"], row["Species"]
        sid = f"{STUDY}__{event}"
        matches = [r for r in s1 if r["SampleID"] == event]
        assert len(matches) == 1 and matches[0]["Species"] == species, event
        specimen = matches[0]
        path = SOURCE / f"{run}.xml"
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={run}&rettype=full&retmode=xml"
        if args.download:
            content = urllib.request.urlopen(url, timeout=60).read()
            tree = ET.fromstring(content)
            public = tree.find(f".//RUN[@accession='{run}']")
            if public is None or public.get("is_public") != "true":
                raise ValueError(f"{run} is not public; no analysis allowed")
            path.write_bytes(content)
        tree = ET.parse(path)
        public = tree.find(f".//RUN[@accession='{run}']")
        assert public is not None and public.get("is_public") == "true", run
        assert tree.find(".//SAMPLE").get("alias") == event
        sample = tree.find(".//SAMPLE")
        assert (
            sample.find("./IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']").text
            == row["BioSample"]
        )
        attributes = {
            e.findtext("TAG"): e.findtext("VALUE")
            for e in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
        }
        project = tree.find(
            ".//STUDY/IDENTIFIERS/EXTERNAL_ID[@namespace='BioProject']"
        ).text
        inputs.append(path)
        paths = []
        for mate in (1, 2):
            pattern = f"maciej_jan_23/Sample_{event}_*/{event}_*_R{mate}_001.fastq.gz"
            if event == "T676-SS3":
                pattern = (
                    f"formalin_ext1/Sample_111936_Barg1_238_339/*_R{mate}_001.fastq.gz"
                )
            candidates = sorted(RAW.glob(pattern))
            assert len(candidates) == 1, (event, mate, candidates)
            raw = candidates[0]
            originals = [
                f
                for f in public.findall("./SRAFiles/SRAFile")
                if f.get("supertype") == "Original"
                and f.get("filename") == f"{event}_R{mate}_001.fastq.gz"
            ]
            assert len(originals) == 1, (run, mate)
            original = originals[0]
            assert raw.stat().st_size == int(original.get("size")), raw
            files.append(
                {
                    "analysis_sample_id": sid,
                    "sra_run": run,
                    "mate": mate,
                    "path": str(raw),
                    "bytes": raw.stat().st_size,
                    "ncbi_original_md5": original.get("md5"),
                    "ncbi_original_filename": original.get("filename"),
                }
            )
            paths.append(str(raw))
        with gzip.open(paths[0], "rt") as handle:
            read_name = handle.readline().strip().lstrip("@").split()[0]
        voucher_s1 = specimen["Catalog"].replace("YPM No. ", "YPM-IZ-")
        voucher_s2 = row["Catalog"]
        notes = []
        voucher = (
            voucher_s2
            if voucher_s1 == voucher_s2 and voucher_s2.startswith("YPM-IZ-")
            else ""
        )
        if voucher_s1 != voucher_s2:
            notes.append(
                f"Voucher unresolved: Table S1 gives {voucher_s1}; Table S2 gives {voucher_s2}. Awaiting author correction; library identity is established by collecting event and public raw-file checksums."
            )
        elif "personal" in voucher_s1:
            notes.append(
                "Specimen in MKM personal collection; no museum voucher number given."
            )
        tissue = attributes.get("tissue", "")
        if not tissue:
            tissue = "nectophore" if event == "T676-SS3" else "unspecified tissue"
            notes.append(
                "Specific tissue not recorded in public metadata; paper describes frozen tissue samples."
                if event != "T676-SS3"
                else "Nectophore specified in paper methods."
            )
        if event == "D725-D8":
            notes.append(
                "Paper methods also mention D725-D9; Tables S1/S2, public SRA and the raw files consistently identify D725-D8."
            )
        metadata.append(
            {
                "analysis_sample_id": sid,
                "library_id": f"{STUDY}:{event}",
                "specimen_id": event,
                "study": STUDY,
                "species_current": species,
                "species_as_published": species,
                "ncbi_species": sample.findtext("./SAMPLE_NAME/SCIENTIFIC_NAME"),
                "collection_id": event,
                "sample_id": "",
                "voucher_id": voucher,
                "voucher_id_s1": voucher_s1,
                "voucher_id_s2": voucher_s2,
                "biosample": row["BioSample"],
                "sra_run": run,
                "bioproject": project,
                "tissue": tissue,
                "preservation": "formalin" if event == "T676-SS3" else "frozen",
                "latitude": f"{float(specimen['Latitude']):.8f}".rstrip("0").rstrip(
                    "."
                ),
                "longitude": f"{float(specimen['Longitude']):.8f}".rstrip("0").rstrip(
                    "."
                ),
                "depth_m": f"{float(specimen['Depth [m]']):.6f}".rstrip("0").rstrip(
                    "."
                ),
                "collection_date": (
                    datetime(1899, 12, 30) + timedelta(days=float(specimen["Date"]))
                )
                .date()
                .isoformat(),
                "location": attributes.get("geo_loc_name", "Pacific Ocean"),
                "ocean_region": "Central Pacific"
                if event == "M120-SS3"
                else "Gulf of California"
                if event == "D725-D8"
                else "NE Pacific",
                "host_route": "none",
                "read_pairs": int(public.get("total_spots")),
                "r1_paths": paths[0],
                "r2_paths": paths[1],
                "sequencing_batches": ":".join(read_name.split(":")[:3]),
                "include_primary": "true",
                "sra_published": public.get("published"),
                "source_doi": DOI,
                "table_s1_row": specimen["source_row"],
                "table_s2_row": row["source_row"],
                "notes": " ".join(notes),
            }
        )
    assert len({r["analysis_sample_id"] for r in metadata}) == 8
    write(OUT / "sample_metadata.tsv", metadata)
    write(OUT / "raw_files.tsv", files)
    write(Path("config/bargmannia_samples.tsv"), metadata)
    # Keep the study manifest compatible with the established cohort workflow.
    inputs.append(Path("manifest.csv"))
    with Path("manifest.csv").open() as handle:
        fields = csv.DictReader(handle).fieldnames
    with Path("manko2026_manifest.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in metadata:
            record = {key: row.get(key, "") for key in fields}
            record.update(
                study_memberships=STUDY,
                original_label=row["collection_id"],
                provenance=f"Manko et al. 2026; doi:{DOI}; Tables S1/S2 and public SRA",
                host_reference="none",
                n_lanes=1,
                locality=row["location"],
                specimen_voucher=row["voucher_id"],
                collection_depth_source="Manko2026 Table S1",
                collection_depth_basis="reported_point",
                depth_source="Table S1",
            )
            writer.writerow(record)
    (OUT / "provenance.json").write_text(
        json.dumps(
            {
                "checked_utc": datetime.now(UTC).isoformat(),
                "source_doi": DOI,
                "selection": "Table S2: Source == This study and raw SRR accession present",
                "inputs_sha256": {str(p): sha(p) for p in inputs},
                "raw_file_identity": "Public NCBI original-upload sizes match local files; MD5 verification required by workflow before processing",
                "species_authority": "Published Table S2, checked against Table S1; original SRA label retained separately",
            },
            indent=2,
        )
        + "\n"
    )
    print(
        json.dumps(
            {
                "libraries": len(metadata),
                "species": dict(Counter(r["species_current"] for r in metadata)),
                "raw_pairs": sum(r["read_pairs"] for r in metadata),
                "capped_pairs": sum(min(r["read_pairs"], 200000000) for r in metadata),
                "tissues": {r["collection_id"]: r["tissue"] for r in metadata},
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
