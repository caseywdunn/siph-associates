"""Explicit review decisions must never silently broaden a metadata correction."""

import importlib
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def review(monkeypatch):
    monkeypatch.syspath_prepend(str(ROOT / "scripts"))
    return importlib.import_module("review_sample_metadata")


def candidate(value="YPM-IZ-35039", **changes):
    return {
        "library_id": "Ahuja2024:CWD16",
        "field": "specimen_voucher",
        "value": value,
        "source_type": "primary_sheet",
        "source": "new.xlsx",
        "source_locator": "Sheet1!C7",
        "source_record": "CWD16",
        "notes": "original note",
        "candidate_file": "sheet_metadata.tsv",
        "candidate_row": "50",
        "source_value": "",
        **changes,
    }


def decision(**changes):
    return {
        "library_id": "Ahuja2024:CWD16",
        "field": "specimen_voucher",
        "source": "new.xlsx",
        "source_locator": "Sheet1!C7",
        "original_value": "YPM-IZ-35039",
        "action": "replace",
        "replacement_value": "YPM:IZ:35039",
        "reason": "Verified punctuation-only voucher normalization",
        **changes,
    }


def test_exact_review_preserves_untargeted_values_and_original_provenance(review):
    rows = [
        candidate(),
        candidate("-74.02", field="longitude", source="old.xlsx"),
        candidate("0-20 m", field="depth_original"),
    ]
    decisions = [
        decision(),
        decision(
            field="longitude",
            source="old.xlsx",
            original_value="-74.02",
            action="exclude",
            replacement_value="",
            reason="Retain verified higher precision record",
        ),
    ]
    result, audit = review.apply_decisions(rows, decisions)
    assert len(result) == 2 and len(audit) == 2
    assert result[0]["value"] == "YPM:IZ:35039"
    assert result[0]["source_value"] == "YPM-IZ-35039"
    assert result[0]["candidate_file"] == "sheet_metadata.tsv"
    assert result[0]["notes"].startswith("original note; review replace:")
    assert result[1] == rows[2], (
        "raw depth differences are not normalized automatically"
    )
    assert audit[1]["source_value"] == "-74.02" and audit[1]["action"] == "exclude"
    assert rows[0]["value"] == "YPM-IZ-35039" and rows[0]["source_value"] == ""
    assert review.apply_decisions(rows, list(reversed(decisions))) == (result, audit)


@pytest.mark.parametrize(
    "changes,match",
    [
        ({"source_locator": "Sheet1!C8"}, "unmatched"),
        ({"original_value": "YPM-IZ-35038"}, "stale"),
        ({"action": "accept"}, "unknown review action"),
        ({"reason": ""}, "requires a reason"),
        ({"replacement_value": ""}, "replacement_value is required"),
    ],
)
def test_invalid_or_stale_decisions_fail(review, changes, match):
    with pytest.raises(ValueError, match=match):
        review.apply_decisions([candidate()], [decision(**changes)])


def test_ambiguous_rows_and_duplicate_decisions_fail(review):
    with pytest.raises(ValueError, match="ambiguous"):
        review.apply_decisions(
            [candidate(), candidate(candidate_file="duplicate.tsv")], [decision()]
        )
    with pytest.raises(ValueError, match="duplicate review"):
        review.apply_decisions([candidate()], [decision(), decision()])


def test_cli_outputs_are_deterministic_and_never_modify_inputs(review, tmp_path):
    candidates, decisions = tmp_path / "candidates.tsv", tmp_path / "decisions.tsv"
    review.write_table(candidates, review.PROVENANCE_COLUMNS, [candidate()], "\t")
    review.write_table(decisions, review.DECISION_COLUMNS, [decision()], "\t")
    before = [path.read_bytes() for path in (candidates, decisions)]
    command = [
        sys.executable,
        str(ROOT / "scripts/review_sample_metadata.py"),
        "--candidates",
        str(candidates),
        "--decisions",
        str(decisions),
        "--output-dir",
    ]
    outputs = [tmp_path / "review1", tmp_path / "review2"]
    for directory in outputs:
        subprocess.run(
            [*command, str(directory)], check=True, capture_output=True, text=True
        )
    for name in ("reviewed_candidates.tsv", "review_decision_audit.tsv"):
        assert (outputs[0] / name).read_bytes() == (outputs[1] / name).read_bytes()
    assert [path.read_bytes() for path in (candidates, decisions)] == before
    repeat = subprocess.run(
        [*command, str(outputs[0])], check=False, capture_output=True, text=True
    )
    assert repeat.returncode != 0
    malformed = tmp_path / "malformed.tsv"
    malformed.write_text("action\nexclude\n")
    with pytest.raises(ValueError, match="decision columns"):
        review.read_decisions(malformed)
