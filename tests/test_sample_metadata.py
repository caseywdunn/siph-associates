"""Collection metadata regressions use tiny fixtures and never open FASTQs."""

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "merge_sample_metadata", ROOT / "scripts/merge_sample_metadata.py"
)
metadata = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(metadata)


def candidate(field, value, library_id="study:one", **provenance):
    return {
        "library_id": library_id,
        "field": field,
        "value": value,
        "source_type": "primary_sheet",
        "source": "sheet.xlsx",
        "source_locator": "Samples!C4",
        "source_record": "one",
        "notes": "original source",
        "candidate_file": "candidates.tsv",
        "candidate_row": "2",
        **provenance,
    }


def test_fill_missing_preserves_immutable_cells_and_provenance():
    rows = [
        {
            "library_id": "study:one",
            "depth_m": "",
            "latitude": "NA",
            "collection_date": "TODO:Table_S7",
            "read_pairs": "123",
            "depth_source": "SRA_spots",
            "sequencing_batches": "instrument:run:flowcell:L001",
            "r1_paths": "/original/R1.fastq.gz",
            "species_current": "Nanomia sp.",
            "include_primary": "true",
            "notes": "keep me",
        }
    ]
    candidates = [
        candidate("depth_m", "0"),
        candidate("latitude", "-31.50"),
        candidate("collection_date", "2024-03-15"),
    ]
    result, accepted, audit = metadata.merge_metadata(rows, candidates)
    assert result[0]["depth_m"] == "0"
    assert result[0]["latitude"] == "-31.5"
    assert result[0]["collection_date"] == "2024-03-15"
    assert {r["status"] for r in audit} == {"filled"}
    assert accepted[0]["source_locator"] == "Samples!C4"
    assert all(r["source"] == "sheet.xlsx" for r in accepted)
    assert metadata.apply_updates(rows, accepted) == result
    for field in rows[0].keys() - {"depth_m", "latitude", "collection_date"}:
        assert result[0][field] == rows[0][field]
    assert rows[0]["depth_m"] == "", "source records must not be mutated"


def test_conflicts_are_held_even_when_source_priority_differs():
    rows = [{"library_id": "study:one", "depth_m": ""}]
    result, accepted, audit = metadata.merge_metadata(
        rows,
        [
            candidate("depth_m", "10"),
            candidate(
                "depth_m", "20", source_type="biosample", source="https://sample/SAM1"
            ),
        ],
    )
    assert result == rows
    assert accepted == []
    assert {r["status"] for r in audit} == {"candidate_conflict"}


def test_agreeing_numeric_formats_keep_all_source_evidence():
    rows = [{"library_id": "study:one", "depth_m": ""}]
    result, accepted, audit = metadata.merge_metadata(
        rows,
        [
            candidate("depth_m", "10.0", source_type="gbif", source="https://gbif/1"),
            candidate("depth_m", "10"),
        ],
    )
    assert result[0]["depth_m"] == "10"
    assert len(accepted) == 2
    assert accepted[0]["source_type"] == "primary_sheet"
    assert {row["value"] for row in accepted} == {"10", "10.0"}
    assert {r["status"] for r in audit} == {"filled"}


def test_existing_values_and_coordinates_cannot_be_overwritten_or_mixed():
    rows = [{"library_id": "study:one", "latitude": "31.000", "longitude": ""}]
    result, accepted, audit = metadata.merge_metadata(
        rows,
        [
            candidate("latitude", "32"),
            candidate("longitude", "170"),
        ],
    )
    assert result == rows
    assert accepted == []
    assert [r["status"] for r in audit] == [
        "preserved_existing_conflict",
        "related_field_conflict",
    ]


def test_ranges_remain_ranges_and_exact_depth_zero_is_retained():
    rows = [{"library_id": "study:one", "depth_m": "0"}]
    result, accepted, _ = metadata.merge_metadata(
        rows,
        [
            candidate("depth_min_m", "0"),
            candidate("depth_max_m", "20"),
            candidate("depth_original", "0-20 m"),
        ],
    )
    assert result[0]["depth_m"] == "0"
    assert result[0]["depth_min_m"] == "0"
    assert result[0]["depth_max_m"] == "20"
    empty_depth = [{"library_id": "study:one", "depth_m": ""}]
    assert metadata.apply_updates(empty_depth, accepted)[0]["depth_m"] == ""


@pytest.mark.parametrize(
    "candidates",
    [
        [candidate("depth_min_m", "30"), candidate("depth_max_m", "20")],
        [
            candidate("depth_m", "25"),
            candidate("depth_min_m", "0"),
            candidate("depth_max_m", "20"),
        ],
    ],
)
def test_inconsistent_depth_fields_are_held_together(candidates):
    rows = [{"library_id": "study:one", "depth_m": ""}]
    result, accepted, audit = metadata.merge_metadata(rows, candidates)
    assert result == rows
    assert accepted == []
    assert {r["status"] for r in audit} == {"related_field_conflict"}


def test_exact_ids_invalid_units_and_immutable_fields_are_rejected():
    rows = [{"library_id": "study:one", "depth_source": "counted", "read_pairs": "12"}]
    result, accepted, audit = metadata.merge_metadata(
        rows,
        [
            candidate("depth_m", "10", library_id="one"),
            candidate("depth_source", "collection sheet"),
            candidate("read_pairs", "100"),
            candidate("latitude", "91"),
            candidate("depth_m", "10-20 m"),
            candidate("collection_date", "2023-02-29"),
        ],
    )
    assert result == rows
    assert accepted == []
    assert [r["status"] for r in audit] == [
        "unknown_library",
        "immutable_field",
        "immutable_field",
        "invalid_candidate",
        "invalid_candidate",
        "invalid_candidate",
    ]


def test_different_original_depth_strings_do_not_block_compatible_measurements():
    rows = [{"library_id": "study:one", "depth_m": ""}]
    result, accepted, audit = metadata.merge_metadata(
        rows,
        [
            candidate("depth_m", "10"),
            candidate("depth_min_m", "0"),
            candidate("depth_max_m", "20"),
            candidate("depth_original", "10 m"),
            candidate("depth_original", "0-20 m", source="second-sheet.xlsx"),
        ],
    )
    assert result[0]["depth_m"] == "10"
    assert result[0]["depth_min_m"] == "0"
    assert result[0]["depth_max_m"] == "20"
    assert "depth_original" not in result[0]
    assert len(accepted) == 3
    assert [r["status"] for r in audit[-2:]] == [
        "candidate_conflict",
        "candidate_conflict",
    ]


def test_raw_coordinate_strings_do_not_block_consistent_numeric_coordinates():
    rows = [
        {
            "library_id": "study:one",
            "latitude": "",
            "longitude": "",
            "lat_long_raw": "37.43 N 72.68 W",
        }
    ]
    result, _, audit = metadata.merge_metadata(
        rows,
        [
            candidate("latitude", "37.43"),
            candidate("longitude", "-72.68"),
            candidate("lat_long_raw", "37.43, -72.68"),
        ],
    )
    assert result[0]["latitude"] == "37.43"
    assert result[0]["longitude"] == "-72.68"
    assert result[0]["lat_long_raw"] == "37.43 N 72.68 W"
    assert [row["status"] for row in audit] == [
        "filled",
        "filled",
        "preserved_existing_conflict",
    ]


def test_different_source_urls_are_combined_without_biological_conflict():
    rows = [{"library_id": "study:one", "depth_m": ""}]
    candidates = [
        candidate("depth_m", "10"),
        candidate("collection_depth_source", "sheet.xlsx#C4"),
        candidate(
            "collection_depth_source",
            "https://sample/SAM1",
            source_type="biosample",
            source="https://sample/SAM1",
        ),
    ]
    result, accepted, audit = metadata.merge_metadata(rows, candidates)
    assert result[0]["depth_m"] == "10"
    assert result[0]["collection_depth_source"] == "sheet.xlsx#C4 | https://sample/SAM1"
    assert {row["status"] for row in audit} == {"filled"}
    source_rows = [row for row in accepted if row["field"] == "collection_depth_source"]
    assert {row["source_value"] for row in source_rows} == {
        "sheet.xlsx#C4",
        "https://sample/SAM1",
    }
    assert metadata.apply_updates(rows, accepted) == result


def test_freeze_captures_nested_metadata_inputs_but_not_historical_freezes(tmp_path):
    spec = importlib.util.spec_from_file_location(
        "freeze_manifest", ROOT / "scripts/freeze_manifest.py"
    )
    freeze = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(freeze)
    relative = [
        "data/sources/original.xlsx",
        "data/sources/metadata_enrichment/sheet_metadata.tsv",
        "data/sources/metadata_enrichment/raw/biosample.json",
        "data/metadata/sample_metadata_updates.tsv",
        "scripts/review_sample_metadata.py",
        "data/metadata/history/old/manifest.freeze.sha256",
    ]
    for name in relative:
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture")
    paths = {
        str(path.relative_to(tmp_path)) for path in freeze.freeze_input_paths(tmp_path)
    }
    assert set(relative[:5]).issubset(paths)
    assert "scripts/merge_sample_metadata.py" in paths
    assert relative[-1] not in paths


def test_rebuild_detects_stale_overlay_instead_of_replacing_source():
    with pytest.raises(ValueError, match="conflicts with source"):
        metadata.apply_updates(
            [{"library_id": "study:one", "depth_m": "40"}], [candidate("depth_m", "10")]
        )
    with pytest.raises(ValueError, match="invalid accepted update"):
        metadata.apply_updates(
            [{"library_id": "study:one"}], [candidate("species_current", "new name")]
        )


def test_cli_writes_only_review_products_and_replays_overlay(tmp_path):
    manifest = tmp_path / "manifest.csv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["library_id", "depth_m", "r1_paths"]
        )
        writer.writeheader()
        writer.writerow(
            {"library_id": "study:one", "depth_m": "", "r1_paths": "/do/not/open.fastq"}
        )
    original_bytes = manifest.read_bytes()
    candidates = tmp_path / "candidates.tsv"
    metadata.write_table(
        candidates,
        metadata.CANDIDATE_COLUMNS,
        [candidate("depth_min_m", "0"), candidate("depth_max_m", "20")],
        "\t",
    )
    output = tmp_path / "review"
    command = [
        sys.executable,
        str(ROOT / "scripts/merge_sample_metadata.py"),
        "--manifest",
        str(manifest),
        "--candidates",
        str(candidates),
        "--updates",
        str(tmp_path / "not-yet-promoted.tsv"),
        "--output-dir",
        str(output),
        "--expected-libraries",
        "1",
    ]
    subprocess.run(command, check=True, capture_output=True, text=True)
    assert manifest.read_bytes() == original_bytes
    result_fields, result = metadata.read_table(output / "manifest.csv", ",")
    accepted = metadata.read_candidates(output / "sample_metadata_updates.tsv")
    assert (
        result_fields
        == ["library_id", "depth_m", "r1_paths"] + metadata.ADDITIONAL_FIELDS
    )
    _, original_rows = metadata.read_table(manifest, ",")
    replayed = metadata.apply_updates(original_rows, accepted)
    assert result == [
        {field: row.get(field, "") for field in result_fields} for row in replayed
    ]
    assert result[0]["depth_m"] == ""
    assert result[0]["depth_min_m"] == "0" and result[0]["depth_max_m"] == "20"
    assert accepted[0]["candidate_file"] == str(candidates)
    assert accepted[0]["source_locator"] == "Samples!C4"
    repeat = subprocess.run(command, check=False, capture_output=True, text=True)
    assert repeat.returncode != 0, "review directories must never be overwritten"
