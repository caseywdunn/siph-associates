#!/usr/bin/env python3
"""Check exclusions, assembly-only changes and preserved focal regional tests."""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

BASE = Path("data/results/coauthor_review/2026-10-22/phase6_analysis/v3")
NEW = Path("data/results/manuscript_revision/2026-10-24/phase6_analysis/v3")
OUT = Path("figures/manuscript_revision")
OUT.mkdir(parents=True, exist_ok=True)
def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
def key(r): return r["sample_id"], r["target_id"]
def present(r, col): return r[col] in {"validated", "high_confidence"}
old = {key(r): r for r in rows(BASE / "grades/bacterial_grades.tsv")}
new = {key(r): r for r in rows(NEW / "grades/bacterial_grades.tsv")}
assert old.keys() == new.keys(), "Complete grade matrix identities changed"
columns = ["grade", "grade_at_5pct", "grade_at_20pct"]
for k in old:
    for col in columns:
        assert present(old[k], col) == present(new[k], col), f"Read-based presence changed: {k}/{col}"
    for col in ["breadth", "reads", "nominated", "mag_from_library"]:
        assert old[k][col] == new[k][col], f"Frozen evidence changed: {k}/{col}"
excluded = {r["sample_id"] for r in rows("config/excluded_libraries.tsv")}
assert len(excluded) == 4 and not any(s.startswith("Church2025") for s in excluded)
for r in new.values():
    assert (r["excluded_library"] == "true") == (r["sample_id"] in excluded)
    assert (r["biological_include"] == "true") == (r["excluded_library"] == r["contaminant_taxon"] == "false")
screen = json.loads(Path("figures/screening/summary.json").read_text())
if "values" in screen: screen = screen["values"]
assert screen["NominatedAccepted"] + screen["OtherLibraryNominatedAccepted"] + screen["MagOnlyPresences"] == screen["ProkaryoticDetections"]
assert screen["NominatedPairs"] == screen["NominatedAccepted"] + screen["NominatedBelow"]
assert screen["NominatedBelow"] == screen["NominatedTrace"] + screen["NominatedNone"]
context = rows("figures/collection_context/library_context.tsv")
assert len(context) == 209 and not excluded & {r["sample_id"] for r in context}
for r in rows("figures/collection_context/associate_detections.tsv"):
    assert r["sample_id"] not in excluded
catalog = rows("figures/manuscript_revision/catalog/bacterial_catalog.manifest.tsv")
cluster = {r["catalog_id"]: r["genome_id"] for r in catalog}
myco = [r for r in new.values() if present(r,"grade") and r["biological_include"] == "true" and "f__Metamycoplasmataceae;" in r["taxonomy"]]
assert len({r["sample_id"] for r in myco}) == 76
assert len({r["target_id"] for r in myco}) == 6
assert not any(cluster[r["target_id"]] == "MAGSP0012" for r in myco)
unchanged = []
empty = []
metadata = {r["sample_id"]:r for r in rows(NEW / "primary/library_flowcells.tsv")}
eligible = {sid for sid,r in metadata.items() if r["host_species"].startswith("Physalia ") and r["inference_eligible"] == "true" and r["ocean_region"]}
assert len(eligible) == 122
for folder,col in [("primary","grade"),("sensitivity/breadth_5pct","grade_at_5pct"),("sensitivity/breadth_20pct","grade_at_20pct")]:
    prior = {r["target_id"]:r for r in rows(BASE / folder / "physalia_region_permutation.tsv")}
    updated = {r["target_id"]:r for r in rows(NEW / folder / "physalia_region_permutation.tsv")}
    for target, row in updated.items():
        if cluster.get(target) not in {"MAGSP0005", "MAGSP0007", "MAGSP0010"}: continue
        for field in ["presences", "deviance_gain", "p_value", "libraries"]:
            assert row[field] == prior[target][field], f"Focal region test changed: {target}/{col}/{field}"
        unchanged.append({"grade_column":col,"target_id":target,"cluster":cluster[target],"p_value":row["p_value"],"old_q_value":prior[target]["q_value"],"q_value":row["q_value"],"focal_statistic_unchanged":True})
    before = {r["sample_id"] for r in old.values() if r["role"] != "decoy" and present(r,col)}
    after = {r["sample_id"] for r in new.values() if r["role"] != "decoy" and present(r,col) and r["biological_include"] == "true"}
    for sid in sorted(eligible):
        empty.append(dict(grade_column=col,sample_id=sid,previously_nonempty=sid in before,retained_nonempty=sid in after,became_empty=sid in before and sid not in after))
def write(name, records):
    with (OUT / name).open("w") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(records[0]),delimiter="\t",lineterminator="\n")
        writer.writeheader();writer.writerows(records)
write("focal_test_comparison.tsv",unchanged)
write("community_profile_changes.tsv",empty)
summary = {"status":"PASS", "grade_pairs_checked":len(new), "presence_calls_unchanged_at_all_three_thresholds":True,
    "assembly_grade_changes":sum(old[k]["grade"] != new[k]["grade"] for k in old),
    "retained_libraries":len(context), "metamycoplasmataceae_libraries":76, "metamycoplasmataceae_detected_clusters":6,
    "focal_regional_statistics_unchanged":len(unchanged),
    "primary_profiles_became_empty":[r["sample_id"] for r in empty if r["grade_column"] == "grade" and r["became_empty"]],
    "limitations":"No new read mapping, dispersion analysis, exclusion of moderate-quality Physalia, or FCS-GX execution."}
(OUT / "validation.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
