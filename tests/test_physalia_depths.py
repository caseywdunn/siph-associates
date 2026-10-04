"""Protect observed depths and distinguish scored surface depths from records."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from merge_sample_metadata import apply_updates, merge_metadata
from score_physalia_depths import depth_candidates


def test_score_only_physalia_without_any_depth_record():
    rows = [
        {
            "library_id": "missing",
            "species_current": "Physalia utriculus",
            "depth_m": "",
        },
        {"library_id": "point", "species_current": "Physalia physalis", "depth_m": "2"},
        {
            "library_id": "interval",
            "species_current": "Physalia sp.",
            "depth_min_m": "0",
            "depth_max_m": "5",
        },
        {
            "library_id": "raw",
            "species_current": "Physalia sp.",
            "depth_original": "unclear depth",
        },
        {"library_id": "other", "species_current": "Nanomia septata", "depth_m": ""},
    ]
    result, accepted, audit = merge_metadata(rows, depth_candidates(rows))
    indexed = {r["library_id"]: r for r in result}
    assert indexed["missing"]["depth_m"] == "0"
    assert indexed["missing"]["collection_depth_basis"] == "curator_assigned_surface"
    assert indexed["point"]["depth_m"] == "2"
    assert indexed["point"]["collection_depth_basis"] == "source_record"
    assert "depth_m" not in indexed["interval"]
    assert "depth_m" not in indexed["raw"]
    assert indexed["other"]["depth_m"] == ""
    assert all(r["status"] == "filled" for r in audit)
    assert apply_updates(rows, accepted) == result
    assert depth_candidates(result) == [], (
        "A second pass must not reclassify assigned depths"
    )
