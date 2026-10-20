"""No read-count exclusion, unique replacement selection, and valid empty outputs."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "completion", ROOT / "scripts/assembly_completion.py"
)
completion = importlib.util.module_from_spec(spec)
spec.loader.exec_module(completion)


def test_no_floor_and_replacements_are_selected_once(tmp_path):
    base = {
        "host_route": "P_physalis",
        "study": "old",
        "r1_paths": "R1",
        "r2_paths": "R2",
        "read_pairs": 100,
    }
    original = [dict(base, sample_id=s) for s in ("kept", "small", "replacement")]
    new = [dict(base, sample_id="", analysis_sample_id="new", study="new")]
    host = [dict(base, sample_id="replacement", host_route="H_rubrum")]
    settings = {
        "minimum_assembly_input_pairs": 0,
        "empty_assembly_is_valid": True,
        "expected_completion_samples": 2,
        "expected_cohort_samples": 4,
        "work_root": "new_results",
    }
    for key, rows in [
        ("original_samples", original),
        ("halistemma_samples", host),
        ("bargmannia_samples", new),
    ]:
        settings[key] = str(tmp_path / f"{key}.tsv")
        completion.write_table(settings[key], rows, list(rows[0]))
    settings["original_membership"] = str(tmp_path / "membership.tsv")
    members = [
        {
            "sample_id": r["sample_id"],
            "assembly_eligible": str(r["sample_id"] != "small").lower(),
        }
        for r in original
    ]
    completion.write_table(settings["original_membership"], members, list(members[0]))
    for name, rows in [("original", original), ("bargmannia", new)]:
        settings[f"{name}_metadata"] = str(tmp_path / f"{name}_metadata.tsv")
        data = [
            {
                "analysis_sample_id": r.get("analysis_sample_id", r["sample_id"]),
                "species_current": "Current species",
            }
            for r in rows
        ]
        completion.write_table(settings[f"{name}_metadata"], data, list(data[0]))
    for name, rows in [
        ("original", original),
        ("halistemma", host),
        ("bargmannia", new),
    ]:
        settings[f"{name}_root"] = str(tmp_path / name)
        for row in rows:
            sample = row.get("analysis_sample_id", row["sample_id"])
            path = (
                Path(settings[f"{name}_root"])
                / "phase3_cohort/inputs"
                / f"{sample}.json"
            )
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(
                json.dumps(
                    dict(
                        row,
                        sample_id=sample,
                        strategy="reference_depleted",
                        retained_pairs=0 if row["sample_id"] == "small" else 88861,
                    )
                )
            )
    attempts, cohort = completion.plan(settings)
    assert {r["sample_id"] for r in attempts} == {"small", "replacement"}
    assert len(cohort) == len({r["sample_id"] for r in cohort}) == 4
    assert all(
        r["assembly_eligible"] == "true" and r["minimum_assembly_input_pairs"] == 0
        for r in cohort
    )
    replacement = next(r for r in cohort if r["sample_id"] == "replacement")
    assert (
        replacement["host_route"] == "H_rubrum"
        and replacement["result_root"] == "new_results"
    )
    settings["minimum_assembly_input_pairs"] = 100000
    with pytest.raises(ValueError, match="no read-count exclusion"):
        completion.plan(settings)


def test_empty_assembly_passes_shared_validation(tmp_path):
    membership = [
        {
            "sample_id": "small",
            "retained_pairs": "88861",
            "minimum_assembly_input_pairs": "0",
            "assembly_eligible": "true",
        }
    ]
    summary = [
        {
            "sample_id": "small",
            "assembly_status": "no_contigs",
            "checkm2_status": "no_bins",
            "raw_bins": "0",
            "assessed_bins": "0",
            "unassessed_bins": "0",
            "binnable_contigs": "0",
            "viral_contigs": "0",
        }
    ]
    for name, rows in [
        ("membership", membership),
        ("eligibility", membership),
        ("summary", summary),
    ]:
        completion.write_table(tmp_path / f"{name}.tsv", rows, list(rows[0]))
    output = tmp_path / "done"
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow/scripts/validate_assemblies.py"),
            "--membership",
            str(tmp_path / "membership.tsv"),
            "--eligibility",
            str(tmp_path / "eligibility.tsv"),
            "--summary",
            str(tmp_path / "summary.tsv"),
            "--output",
            str(output),
        ],
        check=True,
        capture_output=True,
    )
    assert (
        "status\tPASS" in output.read_text() and "no_contigs\t1" in output.read_text()
    )
