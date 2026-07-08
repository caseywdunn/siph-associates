#!/usr/bin/env python3
"""Build the unified specimen/library manifest for siph_associates.

One row per sequencing library (the analysis unit); `specimen_id` (YPM voucher
where available) groups libraries of the same physical animal, and
`also_in_studies` flags specimens shared across studies. Pools the three source
studies symmetrically:

  Church et al. 2025   151 Physalia          -> data/sources/church_samples_metadata.tsv (+ Table S7)
  Ahuja et al. 2024    32-species skim        -> data/sources/ds1_ahuja2024.tsv
  Ahuja et al. 2026    Nanomia population     -> data/sources/ds2_ahuja2026.tsv (+ routing)

Cells that still need external gathering are written as `TODO:<what>` so the
gaps are explicit. Output: manifest.csv (repo root).

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
    """Decimal lat/long from clean sources; DMS-style Ahuja-2024 strings left raw."""
    s = (s or "").strip()
    if "," in s:  # ds2 / Church "lat, long"
        a, b = [x.strip() for x in s.split(",", 1)]
        try:
            return f"{float(a):.5f}", f"{float(b):.5f}"
        except ValueError:
            return "", ""
    return "", ""  # "39.75 N 70.87 W" -> normalize later (kept in lat_long_raw)


rows = []

# ---- Church et al. 2025 (151 Physalia) ----
for r in read_tsv(os.path.join(SRC, "church_samples_metadata.tsv")):
    sid = r["sample"].strip()
    rows.append(dict(
        library_id=f"Church2025:{sid}", specimen_id=sid, study="Church2025",
        original_label=sid, also_in_studies="", provenance="Church et al. 2025",
        species_current=r.get("species", ""), species_as_published=r.get("species", ""),
        host_reference="P_physalis",
        collection_id=sid, ocean_region=r.get("ocean_region", ""), locality=r.get("location", ""),
        latitude=r.get("latitude", ""), longitude=r.get("longitude", ""), lat_long_raw="",
        collection_date="TODO:Table_S7", depth_m="", sra_run="TODO:Church_BioProject",
        bioproject="TODO:Church2025", raw_path_mccleary="TODO:sc2962/config.yaml",
        n_lanes="", notes="cluster=" + r.get("cluster", ""),
    ))

# ---- Nanomia routing -> host_reference for Ahuja 2026 ----
route = {r["sample_id"]: r["decision"] for r in read_tsv(os.path.join(SRC, "nanomia_routing.tsv"))} \
    if os.path.exists(os.path.join(SRC, "nanomia_routing.tsv")) else {}

# ---- Ahuja et al. 2024 (32-species skim) ----
for r in read_tsv(os.path.join(SRC, "ds1_ahuja2024.tsv")):
    sid = r["sample_id"]; lat, lon = split_ll(r.get("lat_lon", ""))
    sp = r.get("species", "")
    href = "N_septata" if ("Nanomia" in sp and "bijuga" not in sp) \
        else "P_physalis" if "Physalia" in sp else "none"
    rows.append(dict(
        library_id=f"Ahuja2024:{sid}", specimen_id=(r.get("voucher") or r.get("collection_id") or sid),
        study="Ahuja2024", original_label=sid,
        also_in_studies=("Ahuja2026" if sid in ("CWD16",) else ""),
        provenance="Ahuja et al. 2024",
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
    sid = r["sample_id"]; lat, lon = split_ll(r.get("lat_lon", ""))
    dec = route.get(sid, "")
    sp = "Nanomia septata" if dec == "R" else "Nanomia sp. (congener; TODO:CO1)"
    href = "N_septata" if dec == "R" else "none"
    rows.append(dict(
        library_id=f"Ahuja2026:{sid}", specimen_id=(r.get("voucher") or sid),
        study="Ahuja2026", original_label=sid,
        also_in_studies=("Ahuja2024" if sid == "NA19" else ""),
        provenance=r.get("source", "Ahuja et al. 2026"),
        species_current=sp, species_as_published="Nanomia (CO1 per Ahuja 2026)",
        host_reference=href, collection_id=r.get("collection_id", ""),
        ocean_region=r.get("ocean", ""), locality=r.get("location", ""),
        latitude=lat, longitude=lon, lat_long_raw=r.get("lat_lon", ""),
        collection_date=r.get("date", ""), depth_m="",
        sra_run=r.get("srr", ""),
        bioproject=("PRJNA925656" if sid == "NA19" else "PRJNA1252167"),
        raw_path_mccleary=r.get("sample_dir", ""), n_lanes=r.get("n_lanes", ""),
        notes=(f"routing properly-paired decision={dec}" if dec else ""),
    ))

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    for r in rows:
        w.writerow({c: r.get(c, "") for c in COLS})

from collections import Counter
print(f"wrote {len(rows)} libraries -> {OUT}")
print("by study:", dict(Counter(r["study"] for r in rows)))
print("by host_reference:", dict(Counter(r["host_reference"] for r in rows)))
todo = sum(1 for r in rows for c in COLS if str(r.get(c, "")).startswith("TODO"))
print(f"TODO cells to fill: {todo}")
