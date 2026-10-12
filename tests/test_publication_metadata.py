"""Publication curation must preserve analysis joins and explicit depth evidence."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "export_sample_metadata", ROOT / "scripts/export_sample_metadata.py"
)
metadata = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(metadata)


def inputs(**values):
    row = {
        "library_id": "study:one",
        "specimen_id": "legacy-specimen",
        "species_current": "Host species",
        "notes": "",
        "locality": "Pacific Ocean",
        "specimen_voucher": "",
        "depth_m": "",
        "depth_min_m": "",
        "depth_max_m": "",
        "collection_depth_basis": "source_record",
        "tissue": "",
        **values,
    }
    samples = [
        {
            "library_id": "study:one",
            "specimen_id": "legacy-specimen",
            "sample_id": "study__one",
        }
    ]
    return [row], samples


def update(field, value, original=""):
    return {
        "library_id": "study:one",
        "field": field,
        "value": value,
        "original_value": original,
        "source": "sheet.tsv",
        "source_locator": "row 2",
    }


def test_biological_ids_do_not_change_workflow_handles_or_input_rows():
    manifest, samples = inputs()
    changes = [
        update("sample_id", "CWD16"),
        update("voucher_id", "YPM:IZ:35039"),
        update("collection_id", "D614-SS12"),
    ]
    result = metadata.curate_metadata(manifest, samples, [], changes)[0]
    assert result["sample_id"] == "CWD16"
    assert result["analysis_sample_id"] == "study__one"
    assert result["library_id"] == "study:one"
    assert result["specimen_id"] == "legacy-specimen"
    assert result["specimen_voucher"] == result["voucher_id"] == "YPM:IZ:35039"
    assert "sample_id" not in manifest[0]
    assert samples[0]["sample_id"] == "study__one"


def test_midpoint_preserves_bounds_and_zero_point_is_not_replaced():
    manifest, samples = inputs(
        depth_min_m="0", depth_max_m="21", depth_original="0-21 m"
    )
    result = metadata.curate_metadata(manifest, samples, [], [])[0]
    assert result["depth_m"] == "10.5"
    assert result["collection_depth_basis"] == "interval_midpoint"
    assert result["depth_min_m"] == "0"
    assert result["depth_max_m"] == "21"
    assert result["depth_original"] == "0-21 m"
    manifest[0]["depth_m"] = "0"
    result = metadata.curate_metadata(manifest, samples, [], [])[0]
    assert result["depth_m"] == "0"
    assert result["collection_depth_basis"] == "source_record"


def test_one_depth_bound_and_unknown_tissue_stay_missing():
    manifest, samples = inputs(depth_min_m="10")
    result = metadata.curate_metadata(manifest, samples, [], [])[0]
    assert result["depth_m"] == ""
    assert result["tissue"] == ""
    assert result["tissue_status"] == "not_recorded"


def test_human_notes_are_separate_from_analysis_annotations():
    manifest, samples = inputs(
        notes="cluster=B2; computational annotation", depth_min_m="0", depth_max_m="20"
    )
    result = metadata.curate_metadata(
        manifest, samples, [], [update("notes", "Collection date needs confirmation.")]
    )[0]
    assert result["analysis_notes"] == manifest[0]["notes"]
    assert "cluster=" not in result["notes"]
    assert result["notes"].startswith("Collection date needs confirmation.")
    assert "midpoint of the reported 0–20 m range" in result["notes"]
    assert "Tissue was not specified" in result["notes"]


@pytest.mark.parametrize(
    "changes",
    [
        [update("library_id", "new:key")],
        [update("specimen_id", "new-specimen")],
        [update("tissue", "tentacle", "stale value")],
        [update("tissue", "tentacle"), update("tissue", "whole body")],
        [{**update("tissue", "tentacle"), "source": ""}],
    ],
)
def test_protected_fields_stale_updates_and_conflicts_are_rejected(changes):
    manifest, samples = inputs()
    with pytest.raises(ValueError):
        metadata.curate_metadata(manifest, samples, [], changes)


@pytest.mark.parametrize(
    "values",
    [
        {"depth_min_m": "20", "depth_max_m": "10"},
        {"depth_min_m": "10", "depth_max_m": "20", "depth_m": "30"},
        {"depth_m": "NaN"},
        {"latitude": "91"},
        {"longitude": "-181"},
    ],
)
def test_invalid_coordinates_and_depths_fail(values):
    manifest, samples = inputs(**values)
    with pytest.raises(ValueError):
        metadata.curate_metadata(manifest, samples, [], [])


def test_duplicate_and_unmatched_workflow_records_fail():
    manifest, samples = inputs()
    for invalid in [
        [],
        samples + samples,
        [{**samples[0], "specimen_id": "different"}],
    ]:
        with pytest.raises(ValueError):
            metadata.curate_metadata(manifest, invalid, [], [])
