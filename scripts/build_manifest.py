#!/usr/bin/env python3
"""Build the unified specimen/library manifest for siph_associates.

One row per sequencing library (the analysis unit); `specimen_id` (YPM voucher
where available) groups libraries of the same physical animal, and
`also_in_studies` flags specimens shared across studies. Pools the three source
studies symmetrically:

  Church et al. 2025   151 Physalia          -> church_samples_metadata.tsv + church_fastq_paths.tsv
  Ahuja et al. 2024    32-species skim        -> ds1_ahuja2024.tsv
  Ahuja et al. 2026    Nanomia population     -> ds2_ahuja2026.tsv + nanomia_routing.tsv

Cells that still need external gathering are written as `TODO:<what>`.
Output: manifest.csv (repo root).

usage: python3 scripts/build_manifest.py
"""
import csv, os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(REPO, "data", "sources")
OUT = os.path.join(REPO, "manifest.csv")

COLS = [
    "library_id", "specimen_id", "study", "original_label", "also_in_studies",
    "provenance", "species_current", "species_as_published", "host_reference",
    "collection_id", "ocean_region", "locality", "latitude", "longitude",
    "lat_long_raw", "collection_date", "depth_m", "sra_run", "bioproject",
    "raw_path_mccleary", "n_lanes", "notes",
]


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def split_ll(s):
    s = (s or "").strip()
    if "," in s:
        a, b = [x.strip() for x in s.split(",", 1)]
        try:
            return f"{float(a):.5f}", f"{float(b):.5f}"
        except ValueError:
            return "", ""
    return "", ""


def norm_voucher(v):
    """YPM vouchers -> NCBI/Darwin-Core colon triplet YPM:IZ:<catalog>; other
    institutions kept as published."""
    v = (v or "").strip()
    if v.upper().startswith("YPM"):
        parts = [p for p in v.replace(":", "-").split("-") if p]
        if len(parts) >= 3:
            return f"YPM:IZ:{parts[-1]}"
    return v


# Ahuja 2026 Nanomia species by locality (paper's CO1 assignments). GoC + France
# are unambiguous; RI and Hawai'i are species-mixed, so the individual call needs
# Ahuja's per-specimen CO1 table.
NANOMIA_SP = {
    "Gulf of California": ("Nanomia bijuga", ""),
    "France": ("Nanomia bijuga", "Villefranche (Mediterranean)"),
    "Rhode Island": ("Nanomia sp.", "RI cluster = N. bijuga (3) + N. cara (1); per-specimen CO1 needed (Ahuja 2026)"),
    "Hawai'i": ("Nanomia sp.", "HI cluster = N. bijuga (2) + N. sp. 1 (1, undescribed); per-specimen CO1 needed (Ahuja 2026)"),
}

church_paths = {r["sample"]: r for r in read_tsv(os.path.join(SRC, "church_fastq_paths.tsv"))} \
    if os.path.exists(os.path.join(SRC, "church_fastq_paths.tsv")) else {}
route = {r["sample_id"]: r["decision"] for r in read_tsv(os.path.join(SRC, "nanomia_routing.tsv"))} \
    if os.path.exists(os.path.join(SRC, "nanomia_routing.tsv")) else {}

rows = []

# ---- Church et al. 2025 (151 Physalia) ----
for r in read_tsv(os.path.join(SRC, "church_samples_metadata.tsv")):
    sid = r["sample"].strip()
    p = church_paths.get(sid, {})
    rows.append(dict(
        library_id=f"Church2025:{sid}", specimen_id=norm_voucher(sid), study="Church2025",
        original_label=sid, also_in_studies="", provenance="Church et al. 2025",
        species_current=r.get("species", ""), species_as_published=r.get("species", ""),
        host_reference="P_physalis",
        collection_id=sid, ocean_region=r.get("ocean_region", ""), locality=r.get("location", ""),
        latitude=r.get("latitude", ""), longitude=r.get("longitude", ""), lat_long_raw="",
        collection_date="TODO:Table_S7", depth_m="", sra_run="TODO:Church_BioProject",
        bioproject="TODO:Church2025", raw_path_mccleary=p.get("dirs", "TODO:sc2962/config.yaml"),
        n_lanes=p.get("n_lanes", ""), notes="cluster=" + r.get("cluster", ""),
    ))

# ---- Ahuja et al. 2024 (32-species skim) ----
for r in read_tsv(os.path.join(SRC, "ds1_ahuja2024.tsv")):
    sid = r["sample_id"]; lat, lon = split_ll(r.get("lat_lon", "")); sp = r.get("species", "")
    href = "N_septata" if ("Nanomia" in sp and "bijuga" not in sp) \
        else "P_physalis" if "Physalia" in sp else "none"
    rows.append(dict(
        library_id=f"Ahuja2024:{sid}", specimen_id=norm_voucher(r.get("voucher") or r.get("collection_id") or sid),
        study="Ahuja2024", original_label=sid,
        also_in_studies=("Ahuja2026" if sid == "CWD16" else ""), provenance="Ahuja et al. 2024",
        species_current=sp, species_as_published=sp, host_reference=href,
        collection_id=r.get("collection_id", ""), ocean_region=r.get("ocean", ""),
        locality=r.get("location", ""), latitude=lat, longitude=lon,
        lat_long_raw=r.get("lat_lon", ""), collection_date=r.get("date", ""), depth_m="",
        sra_run="TODO:Ahuja2024_BioProject", bioproject="TODO:Ahuja2024",
        raw_path_mccleary=r.get("sample_dir", ""), n_lanes=r.get("n_lanes", ""),
        notes=("Physalia utriculus (Guam); maps to P. physalis ref" if sid == "NA22" else ""),
    ))

# ---- Ahuja et al. 2026 (Nanomia population) ----
for r in read_tsv(os.path.join(SRC, "ds2_ahuja2026.tsv")):
    sid = r["sample_id"]; lat, lon = split_ll(r.get("lat_lon", "")); dec = route.get(sid, "")
    loc = r.get("location", "")
    if dec == "R":
        sp, note = "Nanomia septata", "routing R (properly-paired to N. septata ref)"
    else:
        sp, note = NANOMIA_SP.get(loc, ("Nanomia sp.", "congener; per Ahuja 2026 CO1"))
    rows.append(dict(
        library_id=f"Ahuja2026:{sid}", specimen_id=norm_voucher(r.get("voucher") or sid),
        study="Ahuja2026", original_label=sid,
        also_in_studies=("Ahuja2024" if sid == "NA19" else ""),
        provenance=r.get("source", "Ahuja et al. 2026"),
        species_current=sp, species_as_published="Nanomia (CO1 per Ahuja 2026)",
        host_reference=("N_septata" if dec == "R" else "none"),
        collection_id=r.get("collection_id", ""), ocean_region=r.get("ocean", ""),
        locality=loc, latitude=lat, longitude=lon, lat_long_raw=r.get("lat_lon", ""),
        collection_date=r.get("date", ""), depth_m="", sra_run=r.get("srr", ""),
        bioproject=("PRJNA925656" if sid == "NA19" else "PRJNA1252167"),
        raw_path_mccleary=r.get("sample_dir", ""), n_lanes=r.get("n_lanes", ""), notes=note,
    ))

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in COLS})

from collections import Counter
print(f"wrote {len(rows)} libraries -> {OUT}")
print("by study:", dict(Counter(r["study"] for r in rows)))
print("Nanomia species_current:", dict(Counter(r["species_current"] for r in rows if "Nanomia" in r["species_current"])))
todo = sum(1 for r in rows for c in COLS if str(r.get(c, "")).startswith("TODO"))
print(f"TODO cells remaining: {todo}")
