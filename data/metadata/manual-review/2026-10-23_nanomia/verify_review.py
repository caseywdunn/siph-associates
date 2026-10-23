import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

root = Path(__file__).resolve().parents[4]
audit = root / "data/metadata/manual-review/2026-10-23_nanomia"
expected = {
    "Ahuja2026:NA28": ("YPM:IZ:110451", "Nanomia bijuga"),
    "Ahuja2026:WS2": ("YPM:IZ:110216", "Nanomia bijuga"),
    "Ahuja2026:WS3": ("YPM:IZ:111734", "Nanomia bijuga"),
    "Ahuja2026:WS4": ("YPM:IZ:111735", "Nanomia cara"),
    "Ahuja2026:NA29": ("YPM:IZ:111733", "Nanomia sp. 1"),
    "Ahuja2026:NA38": ("YPM:IZ:111591", "Nanomia bijuga"),
}


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


old = {r["library_id"]: r for r in read(audit / "before_sample_metadata.tsv")}
new = {r["library_id"]: r for r in read(root / "figures/cohort/sample_metadata.tsv")}
assert len(new) == 213 and old.keys() == new.keys()
changes = {}
for key, row in new.items():
    diff = {field for field in row if row[field] != old[key][field]}
    if diff:
        assert key in expected and diff == {"species_current", "analysis_notes"}
        changes[key] = sorted(diff)
assert changes.keys() == expected.keys()
for key, (voucher, name) in expected.items():
    assert (new[key]["voucher_id"], new[key]["species_current"]) == (voucher, name)
counts = dict(
    Counter(
        r["species_current"]
        for r in new.values()
        if r["species_current"].startswith("Nanomia")
    )
)
assert counts == {
    "Nanomia bijuga": 10,
    "Nanomia septata": 10,
    "Nanomia cara": 1,
    "Nanomia sp. 1": 1,
}
frozen = json.loads((audit / "before.json").read_text())["sha256"]
unchanged = {}
for path in [
    "manifest.csv",
    "config/samples.tsv",
    "config/post_expansion_samples.tsv",
    "config/post_expansion_manifest.csv",
    "data/metadata/sample_metadata_updates.tsv",
    "figures/cohort/bacterial_detections.tsv",
    "figures/cohort/viral_detections.tsv",
    "figures/cohort/eukaryote_evidence.tsv",
    "figures/coauthor_metadata/community_tests.tsv",
    "figures/coauthor_metadata/dispersion_tests.tsv",
    "figures/cutibacterium_sensitivity/community_tests.tsv",
]:
    current = digest(root / path)
    assert current == frozen[path], path
    unchanged[path] = current
context = read(root / "figures/collection_context/library_context.tsv")
positive = sorted(
    r["library_id"]
    for r in context
    if r["species_current"] == "Nanomia bijuga" and int(r["bacterial_detected"])
)
assert positive == [
    "Ahuja2026:NA31",
    "Ahuja2026:NA33",
    "Ahuja2026:NA37",
    "Ahuja2026:NA38",
]
summary = json.loads((root / "figures/dalila_revision/summary.json").read_text())
assert (
    summary["NanomiaBijugaN"] == 10 and summary["NanomiaBijugaBacterialPositive"] == 4
)
assert summary["NamedSpecies"] == 37 and summary["NamedSpecimens"] == 207
print(
    json.dumps(
        {
            "status": "PASS",
            "date": "2026-10-23",
            "total_libraries": 213,
            "metadata_changes": changes,
            "nanomia_counts": counts,
            "nanomia_bijuga_bacterial_positives": positive,
            "unchanged_sha256": unchanged,
        },
        indent=2,
    )
)
