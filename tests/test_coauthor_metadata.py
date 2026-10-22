"""Protect biological joins and the size of blocked permutation spaces."""

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
SPEC = importlib.util.spec_from_file_location(
    "coauthor_metadata", SCRIPTS / "analyze_coauthor_metadata.py"
)
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def test_tissue_label_allocations_count_distinct_labels_not_read_permutations():
    rows = [
        {"flowcells": "a", "tissue": label}
        for label in ["whole", "whole", "tentacle", "tentacle"]
    ]
    rows += [{"flowcells": "b", "tissue": label} for label in ["whole", "tentacle"]]
    rows += [{"flowcells": "c", "tissue": "tentacle"}] * 3
    assert ANALYSIS.label_allocations(rows, "tissue") == 12


def test_mapping_denominator_requires_same_cohort_even_when_row_counts_match():
    with pytest.raises(ValueError, match="one-to-one"):
        ANALYSIS.join_mapping(
            [{"sample_id": "a"}], [{"sample_id": "b", "input_pairs": "2"}]
        )


def test_duplicate_mapping_denominator_is_rejected():
    rows = [{"sample_id": "a", "input_pairs": "2"}] * 2
    with pytest.raises(ValueError, match="one-to-one"):
        ANALYSIS.join_mapping([{"sample_id": "a"}], rows)


def test_fish_marker_denominator_and_pairs_do_not_replace_genome_mapping_effort():
    rows = [{"sample_id": "a", "input_pairs": 40}]
    fish = [
        {
            "analysis_sample_id": "a",
            "pairs_screened": "90",
            "long_close_fish_match": "2",
        }
    ]
    result = ANALYSIS.join_fish(rows, fish)[0]
    assert result["input_pairs"] == 40
    assert result["fish_12s_pairs_screened"] == 90
    assert result["fish_12s_matching_pairs"] == 2
    assert result["fish_12s_detected"] == 1


def test_depth_design_does_not_treat_california_washington_as_independent_depth():
    rows = [
        {
            "species_current": "Nanomia septata",
            "single_flowcell_eligible": "true",
            "depth_m": str(depth),
            "depth_group": group,
            "locality": locality,
        }
        for depth, group, locality in [
            (10, "shallow", "Washington"),
            (296, "deep", "California"),
        ]
    ]
    result = ANALYSIS.depth_design(rows)[0]
    assert result["status"] == "shallow_deep_contrast_confounded_with_locality"
    assert result["localities_spanning_depth_groups"] == 0
