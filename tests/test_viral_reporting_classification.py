"""Reporting interpretation must not alter mapping evidence or stable target IDs."""

import importlib.util
from pathlib import Path

import pytest

SOURCE = Path(__file__).parents[1] / "scripts/prepare_coauthor_reporting.py"
SPEC = importlib.util.spec_from_file_location("prepare_coauthor_reporting", SOURCE)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_classification_changes_only_interpretation():
    overrides = {
        "VOTU1": {
            "original_class": "associate",
            "reporting_class": "endogenous_candidate",
            "expected_phylum": "Negarnaviricota",
        }
    }
    row = {
        "target_id": "VOTU1",
        "catalog_class": "associate",
        "grade": "high_confidence",
        "reads": "1000",
        "breadth": "0.99",
    }
    result = MODULE.apply_classes([row], overrides, "target_id")
    assert result == [{**row, "catalog_class": "endogenous_candidate"}]
    assert row["catalog_class"] == "associate"


def test_classification_checks_target_and_taxonomy():
    overrides = {
        "VOTU1": {
            "original_class": "associate",
            "reporting_class": "endogenous_candidate",
            "expected_phylum": "Negarnaviricota",
        }
    }
    with pytest.raises(ValueError, match="Missing override"):
        MODULE.apply_classes([], overrides, "target_id")
    with pytest.raises(ValueError, match="Unexpected taxonomy"):
        MODULE.apply_classes(
            [
                {
                    "target_id": "VOTU1",
                    "catalog_class": "associate",
                    "phylum": "Uroviricota",
                }
            ],
            overrides,
            "target_id",
        )
