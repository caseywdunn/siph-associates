#!/usr/bin/env python3
"""Prepare the October 2026 publication metadata additions from retained sources.

Writes a proposal and source-conflict list only. No network or sequence inputs
are accessed. The reviewed proposal is stored separately from analysis metadata.
"""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

from export_sample_metadata import ROOT, missing, read_table, write_table

BASE = Path("data/sources/metadata_enrichment/2026-10-12-publication")
FIELDS = [
    "library_id",
    "field",
    "original_value",
    "value",
    "source",
    "source_locator",
    "notes",
]
CHURCH_COMMIT = "679ac669b4c44f64384373ca7adf15954f055231"
CHURCH_URL = f"https://github.com/shchurch/Physalia_population_genomics/blob/{CHURCH_COMMIT}/data/"


def key(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def voucher(value: str) -> str:
    """Normalize known museum codes only; collection events are not vouchers."""
    match = re.fullmatch(r"YPM[ :_-]*IZ[ :_-]*0*(\d+)", value.strip(), re.I)
    if match:
        return f"YPM:IZ:{int(match[1])}"
    if re.fullmatch(r"(?:FM|NIWA|TMAG|WAM)[ -](?:[A-Z])?\d+", value.strip()):
        return re.sub(r"[ -]+", "-", value.strip())
    return ""


def church_id(value: str) -> str:
    """Use the study's documented museum transfers, not approximate matching."""
    aliases = {"YPM-IZ-111211": "NIWA-173304", "YPM-IZ-111212": "NIWA-173305"}
    value = aliases.get(value, value)
    if re.fullmatch(r"KM563[2-6]", value):
        value = "TMAG-K" + value[2:]
    return key(value)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = {r["library_id"]: r for r in read_table(ROOT / "manifest.csv")}
    candidates = defaultdict(list)
    conflicts = []

    def add(
        library: str, field: str, value: str, source: str, locator: str, notes: str = ""
    ) -> None:
        if library not in manifest:
            raise ValueError(f"Unknown library: {library}")
        candidates[(library, field)].append(
            {
                "library_id": library,
                "field": field,
                "value": value,
                "source": source,
                "source_locator": locator,
                "notes": notes,
            }
        )

    church = {
        church_id(r["original_label"]): r["library_id"]
        for r in manifest.values()
        if r["study"] == "Church2025"
    }
    assert len(church) == 151
    # The paper's specimen table includes institution-qualified vouchers, but
    # the two OA records are collaborator labels without museum attribution.
    for number, row in enumerate(read_table(ROOT / BASE / "church_table_s7.tsv"), 2):
        library = church[church_id(row["ID"])]
        value = voucher(row["ID"])
        if value:
            add(
                library,
                "voucher_id",
                value,
                "data/sources/Church_2025_Table_S7.xlsx",
                f"Table S7!A{number}",
                "Museum identifier from the specimen table; institution prefix retained.",
            )

    # Exact source-row matches were reviewed in the previous enrichment and
    # include the NA22, NA19, and CWD16 aliases across studies.
    matches = {
        (r["source"], r["sheet"], int(r["row"])): r["library_id"]
        for r in read_table(
            ROOT / "data/sources/metadata_enrichment/sheet_row_matches.tsv"
        )
        if r["status"] == "matched"
    }
    raw = json.loads(
        (ROOT / "data/sources/metadata_enrichment/sheet_raw_rows.json").read_text()
    )
    for row in raw:
        library = matches.get((row["source"], row["sheet"], row["row"]))
        if not library:
            continue
        values = row["values"]
        mapping = (
            {
                "sample_id": "sample_name",
                "collection_id": "sample_title",
                "voucher_id": "specimen_voucher",
            }
            if "sample_name" in values
            else {
                "sample_id": "sample_ID",
                "collection_id": "collection_ID",
                "voucher_id": "YPM ID",
            }
        )
        for field, column in mapping.items():
            value = values.get(column, "")
            if missing(value):
                continue
            if field == "voucher_id":
                value = voucher(value)
                if not value:
                    raise ValueError(f"Unrecognized source voucher: {library}")
            add(
                library,
                field,
                value,
                row["source"],
                f"{row['sheet']}!{row['cells'][column]}",
                "Exact reviewed DNA sample match; RNA and out-of-cohort rows excluded.",
            )

    for number, row in enumerate(
        read_table(ROOT / BASE / "physalia_sample_ids.tsv"), 2
    ):
        library = church.get(church_id(row["ID"]))
        if not library:
            continue
        value = row["collection_ID"].strip()
        if missing(value):
            value = ""
        add(
            library,
            "collection_id",
            value,
            CHURCH_URL + "sample_ids.tsv",
            f"row {number};ID={row['ID']};collection_ID",
            "Field event label; a missing event is not replaced by the museum voucher.",
        )
        site = row["site"].strip()
        if not missing(site):
            add(
                library,
                "collection_site",
                site,
                CHURCH_URL + "sample_ids.tsv",
                f"row {number};ID={row['ID']};site",
            )
            location = row["location"].strip()
            value = location if location == site else f"{location}: {site}"
            add(
                library,
                "location",
                value,
                CHURCH_URL + "sample_ids.tsv",
                f"row {number};ID={row['ID']};location,site",
                "Broad source location plus named collection site; frozen locality retained separately.",
            )

    for number, row in enumerate(read_table(ROOT / BASE / "physalia_metadata.tsv"), 2):
        library = church.get(church_id(row["ID"]))
        if not library:
            continue
        sample, tube = (row.get(k) or "" for k in ("YCGA_ID", "diss_tube_ID"))
        if tube and not missing(tube):
            add(
                library,
                "tissue_tube_id",
                tube,
                CHURCH_URL + "metadata.tsv",
                f"row {number};ID={row['ID']};diss_tube_ID",
            )
        value = sample or tube
        if value and not missing(value):
            column = "YCGA_ID" if sample else "diss_tube_ID"
            add(
                library,
                "sample_id",
                value,
                CHURCH_URL + "metadata.tsv",
                f"row {number};ID={row['ID']};{column}",
                "Source sample label; tissue-tube label used when no sequencing-center label is given. Repeated sequencing records must agree.",
            )

    # Fresh public SRA records include 18 samples absent from the old run table.
    # Match exact normalized museum IDs from explicit sample attributes, then
    # require WGS and consistency with every already-recorded accession.
    sra_summary = []
    for path in sorted((ROOT / BASE).glob("physalia_sra_*.xml")):
        for package in ET.parse(path).getroot().findall("EXPERIMENT_PACKAGE"):
            sample = package.find("SAMPLE")
            experiment = package.find("EXPERIMENT")
            if (
                experiment.findtext("./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY")
                != "WGS"
            ):
                continue
            attrs = {
                a.findtext("TAG"): a.findtext("VALUE", "")
                for a in sample.findall("./SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
            }
            labels = [
                sample.get("alias", ""),
                attrs.get("isolate", ""),
                attrs.get("specimen_voucher", ""),
            ]
            libraries = {
                church[church_id(label)]
                for label in labels
                if church_id(label) in church
            }
            if not libraries:
                continue
            if len(libraries) != 1:
                raise ValueError(f"Conflicting SRA specimen matches: {libraries}")
            library = libraries.pop()
            biosamples = {
                e.text
                for e in sample.findall("./IDENTIFIERS/EXTERNAL_ID")
                if e.get("namespace") == "BioSample"
            }
            if len(biosamples) != 1:
                raise ValueError(f"Ambiguous SRA BioSample: {library}")
            biosample = biosamples.pop()
            runs = sorted(r.get("accession") for r in package.findall("./RUN_SET/RUN"))
            run = ";".join(runs)
            source = f"https://www.ncbi.nlm.nih.gov/sra/{experiment.get('accession')}"
            for field, value in (
                ("biosample", biosample),
                ("sra_run", run),
                ("tissue", attrs.get("tissue", "")),
            ):
                if missing(value):
                    continue
                old = manifest[library].get(field, "")
                if not missing(old) and old != value:
                    # Additional experiments for the same specimen are not
                    # automatically identified with this sequenced library.
                    conflicts.append(
                        {
                            "library_id": library,
                            "field": field,
                            "value": value,
                            "source": source,
                            "source_locator": str(path.relative_to(ROOT)),
                            "notes": f"Preserved existing value: {old}",
                        }
                    )
                    continue
                add(
                    library,
                    field,
                    value,
                    source,
                    f"{path.relative_to(ROOT)}#EXPERIMENT={experiment.get('accession')};SAMPLE={sample.get('accession')};{field}",
                    "Exact specimen identifier match in WGS SRA package.",
                )
            if missing(manifest[library].get("country", "")):
                place = attrs.get("geo_loc_name", "").split(":")[0]
                if place and not place.endswith(("Ocean", "Sea")):
                    add(
                        library,
                        "country",
                        place,
                        source,
                        f"SAMPLE={sample.get('accession')};geo_loc_name",
                        "Top-level geographic jurisdiction as submitted.",
                    )
            if missing(manifest[library]["sra_run"]):
                note = manifest[library]["notes"].replace(
                    "; submitted to PRJNA1092115, not yet released by NCBI", ""
                )
                add(
                    library,
                    "analysis_notes",
                    note,
                    source,
                    f"EXPERIMENT={experiment.get('accession')}",
                    "Remove obsolete unreleased-accession note in publication export only.",
                )
            sra_summary.append(
                {
                    "library_id": library,
                    "biosample": biosample,
                    "sra_run": run,
                    "tissue": attrs.get("tissue", ""),
                    "source": source,
                }
            )

    # Carry previously documented specimen-level discrepancies into readable
    # notes so the person reviewing the table need not consult the whole audit.
    review = "data/metadata/enrichment/2026-10-04/README.md"
    notes = {
        "Ahuja2026:NA19": "Depth left blank because sources disagree (296 versus 297 m).",
        "Ahuja2026:WS5": "Depth left blank because the meaning and units of the source value 2854.16 have not been confirmed.",
        "Ahuja2026:WS6": "Depth left blank because the meaning and units of the source value 2563.33 have not been confirmed.",
        "Church2025:YPM-IZ-104465": "Longitude retained as east; the supplemental sheet and BioSample report west despite a Guam locality.",
    }
    for sample in ("NA34", "NA35", "NA36"):
        notes[f"Ahuja2026:{sample}"] = (
            "Collection date retained as 2024-06-16 from the study sheet; BioSample lists 2011-06-16."
        )
    for library, note in notes.items():
        add(
            library,
            "notes",
            note,
            review,
            "Unresolved issues",
            "Human-facing explanation of previously reviewed source disagreement.",
        )

    accepted = []
    for (library, field), group in sorted(candidates.items()):
        values = {r["value"] for r in group}
        # One source spells this known event with an underscore; preserve the
        # primary-sheet spelling and retain the alias as an explicit decision.
        if (
            library == "Church2025:YPM-IZ-104465"
            and field == "collection_id"
            and values == {"Guam-104465", "Guam_104465"}
        ):
            for row in group:
                if row["value"] == "Guam_104465":
                    row["notes"] += (
                        " Source spelling Guam_104465 normalized to primary-sheet Guam-104465 for this exact voucher."
                    )
                    row["value"] = "Guam-104465"
            values = {"Guam-104465"}
        if len(values) > 1:
            conflicts.extend(
                {
                    **r,
                    "notes": r["notes"]
                    + " Conflicting source labels; held for review.",
                }
                for r in group
            )
            continue
        original = manifest[library].get(field, "")
        if field == "voucher_id":
            original = manifest[library].get("specimen_voucher", "")
        elif field == "location":
            original = manifest[library]["locality"]
        elif field == "notes":
            original = ""
        elif field == "analysis_notes":
            original = manifest[library]["notes"]
        accepted.extend({**row, "original_value": original} for row in group)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_table(args.output_dir / "publication_metadata_updates.tsv", accepted, FIELDS)
    write_table(
        args.output_dir / "source_conflicts.tsv",
        conflicts,
        [f for f in FIELDS if f != "original_value"],
    )
    write_table(
        args.output_dir / "sra_specimen_matches.tsv",
        sra_summary,
        ["library_id", "biosample", "sra_run", "tissue", "source"],
    )
    print(
        f"Prepared {len(accepted)} evidence rows; {len(conflicts)} conflict evidence rows"
    )


if __name__ == "__main__":
    main()
