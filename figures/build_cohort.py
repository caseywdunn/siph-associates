#!/usr/bin/env python3
"""Export manuscript figures and source tables from accepted cohort outputs.

Run from the repository root with the existing Snakemake Python environment.
No sequence analysis is performed. All 205 libraries contribute to these
descriptive exports, irrespective of eligibility for flowcell inference.
"""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/cohort"
OUT.mkdir(parents=True, exist_ok=True)
SOURCES = []
plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42,
                     "font.family": "DejaVu Sans", "axes.spines.top": False,
                     "axes.spines.right": False})
VALID = {"validated", "high_confidence"}


def read(relative, delimiter="\t"):
    path = ROOT / relative
    SOURCES.append(path)
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def write(name, rows, fields=None):
    if fields is None:
        fields = list(rows[0])
    with (OUT / name).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


manifest = read("manifest.csv", ",")
meta = {r["library_id"].replace(":", "__", 1): r for r in manifest}
assert len(meta) == 205, "Review figure scope if cohort membership changes"
bac = read("data/results/phase6_analysis/v2/grades/bacterial_grades.tsv")
vir = read("data/results/phase6_analysis/v2/grades/viral_grades.tsv")
euk = read("data/results/eukaryote_gate/eukaryote_grades.tsv")
assert all("count_as_detection" in r for r in euk), "Rebuild corrected eukaryote grades first"
mags = read("data/results/phase3_catalog/mags/mag_catalog.tsv")
catalog = read("data/results/phase4_catalog/v2/bacterial_catalog.manifest.tsv")
assembly = read("config/phase3_assembly.tsv")
flow = read("data/results/phase6_analysis/v2/primary/library_flowcells.tsv")
bp = [r for r in bac if r["grade"] in VALID and r["role"] != "decoy"]
vp = [r for r in vir if r["grade"] in VALID and r["catalog_class"] == "associate"]
ep = [r for r in euk if r["count_as_detection"] == "true"]
assert all(r["grade"] in VALID for r in ep)
host_counts = Counter(r["species_current"] for r in manifest)
studies = ["Church2025", "Ahuja2024", "Ahuja2026"]
colors = {"Church2025": "#0072B2", "Ahuja2024": "#E69F00", "Ahuja2026": "#009E73"}
hosts = sorted(host_counts, key=lambda h: (-host_counts[h], h))
host_rows = [{"host": h, "libraries": host_counts[h],
              **{s: sum(r["species_current"] == h and r["study"] == s for r in manifest)
                 for s in studies}} for h in hosts]
write("host_sampling.tsv", host_rows)
public_fields = ["library_id", "specimen_id", "study", "study_memberships", "species_current",
                 "species_as_published", "ocean_region", "locality", "latitude", "longitude",
                 "collection_date", "depth_m", "depth_min_m", "depth_max_m", "depth_original",
                 "collection_depth_source", "collection_depth_basis", "biosample", "collection_method", "tissue",
                 "life_stage", "country", "specimen_voucher", "sra_run", "bioproject", "sequencing_batches",
                 "read_pairs", "include_primary", "notes"]
write("sample_metadata.tsv", manifest, public_fields)
write("library_flowcells.tsv", flow)
write("bacterial_detections.tsv", bp)
write("viral_detections.tsv", vp, list(vir[0]))
write("eukaryote_evidence.tsv", euk)

# Sampling counts use every library, never the restricted statistical subset.
fig, ax = plt.subplots(figsize=(8.2, 9.4))
left = np.zeros(len(hosts))
for study in studies:
    values = np.array([r[study] for r in host_rows])
    ax.barh(range(len(hosts)), values, left=left, color=colors[study], label=study, height=.75)
    left += values
for i, n in enumerate(left):
    ax.text(n + .55, i, str(int(n)), va="center", fontsize=8)
ax.set_yticks(range(len(hosts)), hosts, fontstyle="italic", fontsize=8)
ax.invert_yaxis()
ax.set_xlim(0, max(left) * 1.1)
ax.set_xlabel("Sequencing libraries (one specimen per library)")
ax.legend(frameon=False, loc="lower right")
ax.set_title("205 libraries across 37 host labels", loc="left", fontweight="bold", pad=12)
fig.tight_layout()
save(fig, "sampling")

# Full-cohort incidence of the seven focal MAG species. Denominators are host
# library counts. A species is counted once per library despite multiple evidence rows.
myco = [r for r in catalog if r["genome_id"] in
        {"MAGSP0005", "MAGSP0007", "MAGSP0010", "MAGSP0029", "MAGSP0031", "MAGSP0011", "MAGSP0012"}]
assert len(myco) == 7
lineages = {"MAGSP0005": "Clade A", "MAGSP0007": "Clade B", "MAGSP0010": "Clade B",
            "MAGSP0029": "Clade B", "MAGSP0031": "Clade B", "MAGSP0011": "DT-68",
            "MAGSP0012": "Mycoplasma_K"}
myco_ids = {r["catalog_id"] for r in myco}
myco_bp = [r for r in bp if r["target_id"] in myco_ids]
present = defaultdict(set)
for r in myco_bp:
    present[(r["target_id"], meta[r["sample_id"]]["species_current"])].add(r["sample_id"])
focal_hosts = [h for h in hosts if h.startswith(("Physalia ", "Nanomia ", "Resomia "))]
counts = np.array([[len(present[(g["catalog_id"], h)]) for g in myco] for h in focal_hosts])
denoms = np.array([host_counts[h] for h in focal_hosts])
fig, ax = plt.subplots(figsize=(9.2, 6.0))
im = ax.imshow(counts / denoms[:, None], vmin=0, vmax=1, cmap="Blues", aspect="auto")
for i, h in enumerate(focal_hosts):
    for j in range(len(myco)):
        value = counts[i, j] / denoms[i]
        ax.text(j, i, f"{counts[i,j]}/{denoms[i]}", ha="center", va="center",
                color="white" if value > .6 else "black", fontsize=8)
ax.set_yticks(range(len(focal_hosts)), focal_hosts, fontstyle="italic")
ax.set_xticks(range(len(myco)), [f"{g['genome_id']}\n{lineages[g['genome_id']]}" for g in myco],
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_title("Mycoplasmatales detection across the full sampled host groups", loc="left", pad=12)
fig.colorbar(im, ax=ax, label="Fraction of libraries with validated presence", shrink=.7)
fig.tight_layout()
save(fig, "mycoplasmatales_incidence")
write("mycoplasmatales_incidence.tsv", [
    {"host": h, "catalog_id": g["catalog_id"], "species_cluster": g["genome_id"],
     "lineage": lineages[g["genome_id"]], "present_libraries": len(present[(g["catalog_id"], h)]),
     "total_libraries": host_counts[h]}
    for h in hosts for g in myco])

# Eukaryote units retain heterogeneous taxonomic ranks. Display detections,
# never claim species richness or add overlapping parent and child units.
units = sorted({r["reporting_unit"] for r in ep})
euk_samples = sorted({r["sample_id"] for r in ep}, key=lambda s: (meta[s]["study"], meta[s]["species_current"], s))
matrix = np.zeros((len(units), len(euk_samples)), dtype=int)
unit_index = {v: i for i, v in enumerate(units)}
sample_index = {v: i for i, v in enumerate(euk_samples)}
grade_code = {"validated": 1, "high_confidence": 2}
for r in ep:
    matrix[unit_index[r["reporting_unit"]], sample_index[r["sample_id"]]] = grade_code[r["grade"]]
unit_labels = [u.split(";")[-1] for u in units]
assert len(unit_labels) == len(set(unit_labels)), "Terminal names collide: show full units"
fig, (ax, bar) = plt.subplots(1, 2, figsize=(12, 8.5), gridspec_kw={"width_ratios": [5, 1]}, sharey=True)
cmap = ListedColormap(["#f3f3f3", "#56B4E9", "#0072B2"])
ax.imshow(matrix, cmap=cmap, norm=BoundaryNorm([-.5,.5,1.5,2.5], 3), aspect="auto")
ax.set_yticks(range(len(units)), unit_labels, fontsize=14)
ax.set_xticks(range(len(euk_samples)), [s.replace("Church2025__", "").replace("Ahuja2024__", "").replace("Ahuja2026__", "") for s in euk_samples], rotation=90, fontsize=13)
ax.set_xlabel("Libraries with validated non-host eukaryotic evidence", fontsize=14)
for j, s in enumerate(euk_samples):
    ax.plot(j, -.9, marker="s", color=colors[meta[s]["study"]], markersize=5.5, clip_on=False)
bar.barh(range(len(units)), (matrix == 1).sum(axis=1), color="#56B4E9")
bar.barh(range(len(units)), (matrix == 2).sum(axis=1), left=(matrix == 1).sum(axis=1), color="#0072B2")
bar.tick_params(axis="y", left=False, labelleft=False)
bar.set_xlabel("Libraries", fontsize=14)
bar.tick_params(axis="x", labelsize=14)
ax.legend(handles=[Patch(color="#56B4E9", label="Read or assembly support"),
                   Patch(color="#0072B2", label="Read and assembly support")]
          + [Patch(color=colors[s], label=s) for s in studies],
          loc="upper left", bbox_to_anchor=(0, 1.20), frameon=False, ncol=3, fontsize=13)
fig.tight_layout()
save(fig, "eukaryote_evidence")

retained = [r for r in mags if r["catalog_status"] == "retained"]
summary = {
    "libraries": len(manifest), "host_labels": len(host_counts),
    "singleton_host_labels": sum(n == 1 for n in host_counts.values()),
    "libraries_by_study": dict(Counter(r["study"] for r in manifest)),
    "assembly_eligible": sum(r["assembly_eligible"] == "true" for r in assembly),
    "retained_mags": len(retained), "mag_species": len({r["species_cluster"] for r in retained}),
    "near_complete_mags": sum(r["near_complete"] == "true" for r in retained),
    "bacterial_detection_rows": len(bp), "bacterial_detected_targets": len({r["target_id"] for r in bp}),
    "bacterial_positive_libraries": len({r["sample_id"] for r in bp}),
    "bacterial_high_confidence_rows": sum(r["grade"] == "high_confidence" for r in bp),
    "bacterial_divergent_match_rows": sum(r["divergent_strain"] == "true" for r in bp),
    "associate_viral_detection_rows": len(vp), "associate_viral_targets": len({r["target_id"] for r in vp}),
    "associate_viral_positive_libraries": len({r["sample_id"] for r in vp}),
    "eukaryote_validated_evidence_rows": sum(r["grade"] in VALID for r in euk),
    "eukaryote_nonredundant_detection_rows": len(ep),
    "eukaryote_positive_libraries": len(euk_samples), "eukaryote_reporting_units": len(units),
    "eukaryote_high_confidence_rows": sum(r["grade"] == "high_confidence" for r in ep),
    "eukaryote_roles": dict(Counter(r["role"] for r in ep)),
    "physalia_mycoplasmatales_positive_libraries": len({r["sample_id"] for r in myco_bp if meta[r["sample_id"]]["species_current"].startswith("Physalia ")}),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
SOURCES.append(Path(__file__).resolve())
write("source_checksums.tsv", [{"path": str(p.relative_to(ROOT)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in SOURCES])
(OUT / "source_revision.txt").write_text(subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip() + "\nWorking-tree corrections are recorded by source_checksums.tsv and review archive.\n")
print(json.dumps(summary, indent=2))
