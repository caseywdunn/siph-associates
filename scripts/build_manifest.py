#!/usr/bin/env python3
"""Build publication manifest, paired-FASTQ table, and provenance table.

The analytical manifest has one row per unique raw library. Study provenance is
many-to-one: Church YPM-IZ-104465 and Ahuja 2024 NA22 are the same library,
whereas CWD16 and NA19 are distinct Ahuja 2024 libraries reused by Ahuja 2026.
"""
from __future__ import annotations

import csv
import glob
import gzip
import os
import re
import site
import sys
from collections import Counter, defaultdict
from pathlib import Path

from merge_sample_metadata import ADDITIONAL_FIELDS, apply_updates, read_candidates

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "data" / "sources"
META = REPO / "data" / "metadata"
OUT = REPO / "manifest.csv"
RAW_OUT = META / "raw_files.tsv"
PROV_OUT = META / "library_provenance.tsv"
STUDY_ORDER = {"Church2025": 0, "Ahuja2024": 1, "Ahuja2026": 2}

COLS = [
    "library_id", "specimen_id", "study", "study_memberships",
    "original_label", "also_in_studies", "provenance", "species_current",
    "species_as_published", "host_reference", "collection_id",
    "ocean_region", "locality", "latitude", "longitude", "lat_long_raw",
    "collection_date", "depth_m", "sra_run", "bioproject",
    "sequencing_batches", "raw_path_mccleary", "r1_paths", "r2_paths",
    "n_lanes", "read_pairs", "depth_source", "include_primary",
    "exclusion_reason", "notes", *ADDITIONAL_FIELDS,
]
RAW_COLS = [
    "library_id", "read_pair_id", "sequencing_batch", "source_directory",
    "r1_path", "r2_path", "r1_bytes", "r2_bytes", "r1_read_length",
    "r2_read_length", "sequencing_batch_source",
]
PROV_COLS = [
    "provenance_id", "library_id", "study", "original_label", "specimen_id",
    "record_role", "sequence_origin", "bioproject", "sra_run",
    "source_record", "notes",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader(
            (line for line in fh if not line.startswith("#")), delimiter="\t"
        ))


def write_table(path: Path, rows: list[dict[str, object]], cols: list[str],
                delimiter: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=cols, delimiter=delimiter, lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, "") for col in cols})


def split_ll(value: str) -> tuple[str, str]:
    value = (value or "").strip()
    if "," in value:
        left, right = [part.strip() for part in value.split(",", 1)]
        try:
            return f"{float(left):.5f}", f"{float(right):.5f}"
        except ValueError:
            pass
    return "", ""


def norm_voucher(value: str) -> str:
    value = (value or "").strip()
    if value.upper().startswith("YPM"):
        parts = [part for part in value.replace(":", "-").split("-") if part]
        if len(parts) >= 3:
            return f"YPM:IZ:{parts[-1]}"
    return value


def norm_key(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]", "", (value or "")).upper()


def ordered_unique(values: list[str], study_order: bool = False) -> list[str]:
    found = {value for value in values if value}
    if study_order:
        return sorted(found, key=lambda value: (STUDY_ORDER.get(value, 99), value))
    return sorted(found)


def mate_path(r1: str) -> str:
    if "_R1_" in r1:
        return r1.replace("_R1_", "_R2_", 1)
    if r1.endswith("_R1.fastq.gz"):
        return r1[:-len("_R1.fastq.gz")] + "_R2.fastq.gz"
    raise ValueError(f"cannot derive R2 name from {r1}")


def pair_id(r1: str) -> str:
    name = Path(r1).name
    name = re.sub(r"_R1_\d+\.fastq\.gz$", "", name)
    return re.sub(r"_R1\.fastq\.gz$", "", name)


def sequencing_batch(r1: str, header: str) -> str:
    """Use the actual Illumina instrument/run/flowcell/lane header fields."""
    token = header.removeprefix("@").split()[0]
    fields = token.split(":")
    if len(fields) >= 4 and all(fields[:4]):
        return f"{fields[0]}:{fields[1]}:{fields[2]}:L{fields[3].zfill(3)}"
    # Kept only as an explicit fallback for non-CASAVA data; validation rejects
    # it in the current production manifest so batch never silently degrades.
    path = Path(r1)
    delivery = path.parent.parent.name
    flowcell = re.match(r"^([A-Z0-9]+)_([0-9]+)_", path.name)
    if flowcell:
        return f"{delivery}:{flowcell.group(1)}:lane{flowcell.group(2)}"
    lane = re.search(r"_L(\d{3})(?:_|\.)", path.name)
    if lane:
        return f"{delivery}:L{lane.group(1)}"
    return f"{delivery}:lane_unknown"


def first_read_metadata(path: str) -> tuple[int, str]:
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n")
        sequence = fh.readline().rstrip("\n")
    if not header.startswith("@") or not sequence:
        raise ValueError(f"not a readable FASTQ record: {path}")
    return len(sequence), header


def scan_illumina_r1(directory: str) -> list[str]:
    """Select original lane files and exclude derived combined_R1 files."""
    pattern = os.path.join(directory, "*_R1_[0-9][0-9][0-9].fastq.gz")
    return sorted(glob.glob(pattern))


NANOMIA_SP = {
    "Gulf of California": ("Nanomia bijuga", ""),
    "France": ("Nanomia bijuga", "Villefranche (Mediterranean)"),
    "Rhode Island": (
        "Nanomia sp.",
        "RI cluster = N. bijuga (3) + N. cara (1); per-specimen CO1 needed",
    ),
    "Hawai'i": (
        "Nanomia sp.",
        "HI cluster = N. bijuga (2) + N. sp. 1 (1); per-specimen CO1 needed",
    ),
}

church_paths = {
    row["sample"]: row for row in read_tsv(SRC / "church_fastq_paths.tsv")
}
route = {
    row["sample_id"]: row["decision"]
    for row in read_tsv(SRC / "nanomia_routing.tsv")
}
RUN_SPOTS: dict[str, int] = {}


def load_sra(patterns: list[str]) -> dict[str, tuple[str, str]]:
    """Load header-less SRA runinfo: run=0, spots=3, project=21, sample=29."""
    result: dict[str, tuple[str, str]] = {}
    for pattern in patterns:
        for path in glob.glob(str(SRC / pattern)):
            with open(path) as fh:
                for row in csv.reader(fh):
                    if len(row) < 30 or not row[0].startswith(("SRR", "ERR", "DRR")):
                        continue
                    try:
                        RUN_SPOTS[row[0]] = int(row[3])
                    except (ValueError, IndexError):
                        pass
                    for key in {norm_key(row[29]), norm_key(row[11])}:
                        if not key:
                            continue
                        runs, project = result.get(key, ("", row[21]))
                        accessions = [value for value in runs.split(";") if value]
                        if row[0] not in accessions:
                            accessions.append(row[0])
                        result[key] = (";".join(accessions), row[21])
    return result


def sra_lookup(table: dict[str, tuple[str, str]], *candidates: str):
    for candidate in candidates:
        hit = table.get(norm_key(candidate))
        if hit:
            return hit
    return None


SRA_CHURCH = load_sra(["church*_runinfo.csv"])
SRA_A2024 = load_sra(["ahuja2024*_runinfo.csv"])
load_sra(["ahuja2026*_runinfo.csv"])
COUNTED = {
    row["library_id"]: int(row["read_pairs"])
    for row in read_tsv(SRC / "counted_read_pairs.tsv")
}

# Church collection dates from Supplementary Table S7.
CHURCH_DATE: dict[str, str] = {}
try:
    sys.path.insert(0, site.getusersitepackages())
    import openpyxl  # type: ignore

    workbook = openpyxl.load_workbook(
        SRC / "Church_2025_Table_S7.xlsx", read_only=True, data_only=True
    )
    sheet = workbook.active
    header = None
    date_index = None
    for values in sheet.iter_rows(values_only=True):
        if header is None:
            header = [str(cell or "") for cell in values]
            date_index = next(
                (index for index, name in enumerate(header) if "date" in name.lower()),
                None,
            )
            continue
        if values and values[0] and date_index is not None and values[date_index]:
            CHURCH_DATE[norm_key(str(values[0]))] = str(values[date_index]).split()[0]
except Exception as error:  # noqa: BLE001
    print(f"Church dates unavailable: {error}", file=sys.stderr)

source_rows: list[dict[str, object]] = []
church_pending: list[str] = []

# Church et al. 2025: all 151 published libraries.
for source in read_tsv(SRC / "church_samples_metadata.tsv"):
    sample = source["sample"].strip()
    hit = sra_lookup(
        SRA_CHURCH, sample, norm_voucher(sample).replace("YPM:IZ:", "YPM-IZ-")
    )
    paths = church_paths.get(sample, {})
    note = "cluster=" + source.get("cluster", "")
    if not hit:
        church_pending.append(sample)
        note += "; submitted to PRJNA1092115, not yet released by NCBI"
    source_rows.append(dict(
        library_id=f"Church2025:{sample}", specimen_id=norm_voucher(sample),
        study="Church2025", original_label=sample,
        provenance="Church et al. 2025",
        species_current=source.get("species", ""),
        species_as_published=source.get("species", ""),
        host_reference="P_physalis", collection_id=sample,
        ocean_region=source.get("ocean_region", ""),
        locality=source.get("location", ""), latitude=source.get("latitude", ""),
        longitude=source.get("longitude", ""), lat_long_raw="",
        collection_date=CHURCH_DATE.get(norm_key(sample), "TODO:Table_S7"),
        depth_m="", sra_run=(hit[0] if hit else ""),
        bioproject=(hit[1] if hit else "PRJNA1092115"),
        raw_path_mccleary=paths.get("dirs", "TODO:sc2962/config.yaml"),
        n_lanes=paths.get("n_lanes", ""), notes=note,
        source_record="church_samples_metadata.tsv",
    ))

# Ahuja et al. 2024 broad skim.
for source in read_tsv(SRC / "ds1_ahuja2024.tsv"):
    sample = source["sample_id"]
    lat, lon = split_ll(source.get("lat_lon", ""))
    species = source.get("species", "")
    host_ref = (
        "N_septata" if ("Nanomia" in species and "bijuga" not in species)
        else "P_physalis" if "Physalia" in species else "none"
    )
    hit = sra_lookup(SRA_A2024, sample, source.get("collection_id", ""))
    source_rows.append(dict(
        library_id=f"Ahuja2024:{sample}",
        specimen_id=norm_voucher(
            source.get("voucher") or source.get("collection_id") or sample
        ),
        study="Ahuja2024", original_label=sample,
        provenance="Ahuja et al. 2024", species_current=species,
        species_as_published=species, host_reference=host_ref,
        collection_id=source.get("collection_id", ""),
        ocean_region=source.get("ocean", ""), locality=source.get("location", ""),
        latitude=lat, longitude=lon, lat_long_raw=source.get("lat_lon", ""),
        collection_date=source.get("date", ""), depth_m="",
        sra_run=(hit[0] if hit else "TODO:SRA"),
        bioproject=(hit[1] if hit else "PRJNA925656"),
        raw_path_mccleary=source.get("sample_dir", ""),
        n_lanes=source.get("n_lanes", ""),
        notes=("Physalia utriculus (Guam); maps to P. physalis ref"
               if sample == "NA22" else ""),
        source_record="ds1_ahuja2024.tsv",
    ))

# Ahuja et al. 2026 Nanomia population set.
for source in read_tsv(SRC / "ds2_ahuja2026.tsv"):
    sample = source["sample_id"]
    lat, lon = split_ll(source.get("lat_lon", ""))
    decision = route.get(sample, "")
    locality = source.get("location", "")
    if decision == "R":
        species, note = (
            "Nanomia septata",
            "routing R (properly-paired to N. septata ref)",
        )
    else:
        species, note = NANOMIA_SP.get(
            locality, ("Nanomia sp.", "congener; per Ahuja 2026 CO1")
        )
    source_rows.append(dict(
        library_id=f"Ahuja2026:{sample}",
        specimen_id=norm_voucher(source.get("voucher") or sample),
        study="Ahuja2026", original_label=sample,
        provenance=source.get("source", "Ahuja et al. 2026"),
        species_current=species,
        species_as_published="Nanomia (CO1 per Ahuja 2026)",
        host_reference=("N_septata" if decision == "R" else "none"),
        collection_id=source.get("collection_id", ""),
        ocean_region=source.get("ocean", ""), locality=locality,
        latitude=lat, longitude=lon, lat_long_raw=source.get("lat_lon", ""),
        collection_date=source.get("date", ""), depth_m="",
        sra_run=source.get("srr", ""),
        bioproject=("PRJNA925656" if sample == "NA19" else "PRJNA1252167"),
        raw_path_mccleary=source.get("sample_dir", ""),
        n_lanes=source.get("n_lanes", ""), notes=note,
        source_record="ds2_ahuja2026.tsv",
    ))

# Exact read-pair counts: deposited SRA spots, otherwise prior full FASTQ counts.
for row in source_rows:
    runs = [value for value in str(row.get("sra_run", "")).split(";") if value]
    spots = [RUN_SPOTS[value] for value in runs if value in RUN_SPOTS]
    if runs and len(spots) == len(runs):
        row["read_pairs"] = sum(spots)
        row["depth_source"] = "SRA_spots"
    elif row["library_id"] in COUNTED:
        row["read_pairs"] = COUNTED[str(row["library_id"])]
        row["depth_source"] = "counted"
    else:
        row["read_pairs"] = ""
        row["depth_source"] = "TODO:count"

# Canonical analytical identity. NA22 is byte-identical to the Church library.
canonical = {"Ahuja2024:NA22": "Church2025:YPM-IZ-104465"}
analysis_rows = [
    dict(row) for row in source_rows if row["library_id"] not in canonical
]
analysis_by_id = {str(row["library_id"]): row for row in analysis_rows}

provenance_rows: list[dict[str, object]] = []
for row in source_rows:
    source_id = str(row["library_id"])
    library_id = canonical.get(source_id, source_id)
    provenance_rows.append(dict(
        provenance_id=source_id, library_id=library_id, study=row["study"],
        original_label=row["original_label"], specimen_id=row["specimen_id"],
        record_role=("duplicate_provenance" if source_id in canonical
                     else "primary_analysis"),
        sequence_origin=row["provenance"], bioproject=row["bioproject"],
        sra_run=row["sra_run"], source_record=row["source_record"], notes=row["notes"],
    ))

# Both are old Ahuja 2024 libraries explicitly reused in Ahuja 2026 Table S8.
cwd16 = analysis_by_id["Ahuja2024:CWD16"]
provenance_rows.append(dict(
    provenance_id="Ahuja2026:CWD16", library_id="Ahuja2024:CWD16",
    study="Ahuja2026", original_label="CWD16", specimen_id=cwd16["specimen_id"],
    record_role="reused_in_study", sequence_origin="Ahuja et al. 2024",
    bioproject="PRJNA925656", sra_run="SRR23143286",
    source_record="Ahuja_2026_journal.pone.0351247.s008.xlsx",
    notes="Distinct from NA19; listed as older data in Ahuja 2026 Table S8",
))
na19 = analysis_by_id["Ahuja2026:NA19"]
provenance_rows.append(dict(
    provenance_id="Ahuja2024:NA19", library_id="Ahuja2026:NA19",
    study="Ahuja2024", original_label="NA19", specimen_id=na19["specimen_id"],
    record_role="sequence_origin", sequence_origin="Ahuja et al. 2024",
    bioproject="PRJNA925656", sra_run="SRR23143273",
    source_record="Ahuja_2024_Supplementary_Table_1.xlsx",
    notes="Same library reused in the Ahuja 2026 Nanomia analysis",
))

provenance_by_library: dict[str, list[dict[str, object]]] = defaultdict(list)
for row in provenance_rows:
    provenance_by_library[str(row["library_id"])].append(row)

for row in analysis_rows:
    library_id = str(row["library_id"])
    records = provenance_by_library[library_id]
    memberships = ordered_unique([str(record["study"]) for record in records], True)
    row["study_memberships"] = ";".join(memberships)
    row["also_in_studies"] = ";".join(
        study for study in memberships if study != row["study"]
    )
    row["provenance"] = ";".join(ordered_unique(
        [str(record["sequence_origin"]) for record in records]
    ))
    row["sra_run"] = ";".join(ordered_unique(
        [str(record["sra_run"]) for record in records]
    ))
    row["bioproject"] = ";".join(ordered_unique(
        [str(record["bioproject"]) for record in records]
    ))
    row["include_primary"] = "true"
    row["exclusion_reason"] = ""
    if library_id == "Church2025:YPM-IZ-104465":
        row["notes"] = str(row["notes"]) + "; also Ahuja2024:NA22"
        row["read_pairs"] = RUN_SPOTS["SRR23143271"]
        row["depth_source"] = "SRA_spots"

# Normalize selected raw file pairs and attach aggregate lists to the manifest.
raw_rows: list[dict[str, object]] = []
for row in analysis_rows:
    library_id = str(row["library_id"])
    if row["study"] == "Church2025":
        path_record = church_paths[str(row["original_label"])]
        r1_paths = [value for value in path_record["R1_paths"].split(";") if value]
    else:
        r1_paths = scan_illumina_r1(str(row["raw_path_mccleary"]))
    if not r1_paths:
        raise ValueError(f"no selected R1 FASTQs for {library_id}")

    library_raw: list[dict[str, object]] = []
    for r1 in sorted(r1_paths):
        r2 = mate_path(r1)
        if not Path(r1).is_file() or not Path(r2).is_file():
            raise FileNotFoundError(f"missing pair for {library_id}: {r1} / {r2}")
        r1_length, r1_header = first_read_metadata(r1)
        r2_length, _ = first_read_metadata(r2)
        raw = dict(
            library_id=library_id, read_pair_id=pair_id(r1),
            sequencing_batch=sequencing_batch(r1, r1_header),
            source_directory=str(Path(r1).parent), r1_path=r1, r2_path=r2,
            r1_bytes=Path(r1).stat().st_size, r2_bytes=Path(r2).stat().st_size,
            r1_read_length=r1_length, r2_read_length=r2_length,
            sequencing_batch_source="FASTQ header: instrument/run/flowcell/lane",
        )
        library_raw.append(raw)
        raw_rows.append(raw)

    row["sequencing_batches"] = ";".join(ordered_unique(
        [str(raw["sequencing_batch"]) for raw in library_raw]
    ))
    row["raw_path_mccleary"] = ";".join(ordered_unique(
        [str(raw["source_directory"]) for raw in library_raw]
    ))
    row["r1_paths"] = ";".join(str(raw["r1_path"]) for raw in library_raw)
    row["r2_paths"] = ";".join(str(raw["r2_path"]) for raw in library_raw)
    row["n_lanes"] = len(library_raw)

analysis_rows.sort(key=lambda row: (
    STUDY_ORDER.get(str(row["study"]), 99), str(row["library_id"])
))
raw_rows.sort(key=lambda row: (str(row["library_id"]), str(row["read_pair_id"])))
provenance_rows.sort(key=lambda row: str(row["provenance_id"]))

# Reviewed collection metadata is an overlay on the deduplicated analytical IDs.
# This does not alter FASTQ identity, sequencing depth_source, or frozen results.
metadata_updates = META / "sample_metadata_updates.tsv"
if metadata_updates.exists():
    analysis_rows = apply_updates(analysis_rows, read_candidates(metadata_updates))

write_table(OUT, analysis_rows, COLS, ",")
write_table(RAW_OUT, raw_rows, RAW_COLS, "\t")
write_table(PROV_OUT, provenance_rows, PROV_COLS, "\t")

print(f"wrote {len(analysis_rows)} unique libraries -> {OUT}")
print("by primary study:", dict(Counter(row["study"] for row in analysis_rows)))
print(f"wrote {len(raw_rows)} paired FASTQ records -> {RAW_OUT}")
print(f"wrote {len(provenance_rows)} provenance records -> {PROV_OUT}")
print(f"Church SRA pending release: {len(church_pending)}")
todo = sum(
    1 for row in analysis_rows for col in COLS
    if str(row.get(col, "")).startswith("TODO")
)
print(f"TODO cells remaining: {todo}")
