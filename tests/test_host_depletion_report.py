"""Check sample joins, pair/read units, and incomplete mapping treatment."""

import csv
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "host_report", ROOT / "scripts/summarize_host_depletion.py"
)
report = importlib.util.module_from_spec(spec)
spec.loader.exec_module(report)


def write_table(path, rows):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@pytest.fixture
def mapping(tmp_path):
    sample = "study__specimen"
    metadata = tmp_path / "metadata.tsv"
    write_table(
        metadata,
        [{"analysis_sample_id": sample, "species_current": "Physalia utriculus"}],
    )
    manifest = tmp_path / "samples.tsv"
    write_table(
        manifest,
        [
            {
                "sample_id": sample,
                "host_route": "P_physalis",
                "include_primary": "true",
            },
            {"sample_id": "other", "host_route": "none", "include_primary": "true"},
        ],
    )
    settings = {
        "metadata": str(metadata),
        "sources": [{"samples": str(manifest), "results": str(tmp_path)}],
        "allowed_pending_samples": [],
    }
    source = report.selected_samples(settings)[sample]
    trim = {"sample_id": sample, "raw_pairs_expected": 100, "pairs_after_fastp": 3}
    metrics = {
        "sample_id": sample,
        "host_route": "P_physalis",
        "strategy": "reference_depleted",
        "input_pairs": 3,
        "retained_pairs": 1,
        "retained_fraction": 1 / 3,
    }
    for key, data in (("trim", trim), ("metrics", metrics)):
        path = Path(source[key])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(data))
    provenance = {k: metrics[k] for k in ("sample_id", "host_route", "strategy")} | {
        "upstream_provenance_sha256": report.sha256(source["trim"])
    }
    path = Path(source["mapping"])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(provenance))
    return settings, sample, source


def test_current_species_post_qc_reads_and_pair_removal(mapping):
    settings, sample, _ = mapping
    rows, pending, _ = report.build_rows(settings)
    # Three input pairs: both mates map in one, one mate in another, neither
    # mate in the third. Removed fraction is 2/3 of pairs, not 3/6 mapped reads.
    assert rows == [
        {
            "species": "Physalia utriculus",
            "sample": sample,
            "reads_ingested": 6,
            "percent_removed": "66.6667",
            "reads_remaining": 2,
        }
    ]
    assert not pending


def test_pending_mapping_is_explicit_and_never_zero(mapping):
    settings, sample, source = mapping
    Path(source["mapping"]).unlink()
    with pytest.raises(FileNotFoundError, match="Missing completed"):
        report.build_rows(settings)
    settings["allowed_pending_samples"] = [sample]
    rows, pending, _ = report.build_rows(settings)
    assert rows[0]["reads_ingested"] == 6
    assert rows[0]["reads_remaining"] == rows[0]["percent_removed"] == "NA"
    assert pending == [sample]
    assert source["mapping"] not in report.report_inputs(settings)
    assert source["mapping"] in report.report_inputs(settings, require_complete=True)
    with pytest.raises(FileNotFoundError):
        report.build_rows(settings, require_complete=True)


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("sample_id", "wrong", "identity/strategy"),
        ("host_route", "H_rubrum", "identity/strategy"),
        ("strategy", "fixed_effort", "identity/strategy"),
        ("input_pairs", 4, "counts disagree"),
        ("retained_pairs", 4, "count bounds"),
        ("retained_fraction", 0.2, "fraction disagrees"),
    ],
)
def test_inconsistent_mapping_is_rejected(mapping, field, value, message):
    settings, _, source = mapping
    path = Path(source["metrics"])
    metrics = json.loads(path.read_text())
    metrics[field] = value
    path.write_text(json.dumps(metrics))
    with pytest.raises(ValueError, match=message):
        report.build_rows(settings)


def test_mapping_must_reference_matching_trim_record(mapping):
    settings, _, source = mapping
    path = Path(source["trim"])
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError, match="provenance mismatch"):
        report.build_rows(settings)


def test_metadata_join_must_be_unique_and_complete(mapping):
    settings, sample, _ = mapping
    path = Path(settings["metadata"])
    row = {"analysis_sample_id": sample, "species_current": "Physalia utriculus"}
    write_table(path, [row, row])
    with pytest.raises(ValueError, match="Duplicate analysis_sample_id"):
        report.build_rows(settings)
    write_table(path, [dict(row, analysis_sample_id="different")])
    with pytest.raises(ValueError, match="Missing species metadata"):
        report.build_rows(settings)


def test_completed_pending_sample_automatically_enters_dependencies(mapping):
    settings, sample, source = mapping
    settings["allowed_pending_samples"] = [sample]
    inputs = report.report_inputs(settings)
    assert source["metrics"] in inputs and source["mapping"] in inputs
    _, pending, _ = report.build_rows(settings)
    assert not pending
