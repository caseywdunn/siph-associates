"""Protect specimen identity and distinct bacterial/human evidence units."""

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "human_cutibacterium",
    Path(__file__).resolve().parents[1] / "scripts/summarize_human_cutibacterium.py",
)
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def fixtures():
    metadata = [
        dict(
            analysis_sample_id=sample,
            library_id=sample,
            study="study",
            species_current="species",
            ocean_region="region",
            locality="place",
            tissue="tissue",
        )
        for sample in ("a", "b", "c")
    ]
    human = [
        dict(
            sample_id=sample,
            selection="equal_per_source_file_prefix",
            raw_selected_pairs="2000",
            qc_pairs="1000",
            stringent_human_compatible_pairs=count,
            distinct_sequence_pairs=count,
            reference_sequences="1",
            occupied_1Mb_windows="1",
        )
        for sample, count in [("a", "10"), ("b", "0"), ("c", "1")]
    ]
    bacteria = [
        dict(
            sample_id="a",
            target_id="BAC00019",
            grade="validated",
            role="reference",
            breadth="0.5",
            reads="40",
            assembly_support_bases="0",
        ),
        dict(
            sample_id="b",
            target_id="BAC00019",
            grade="not_detected",
            role="reference",
            breadth="0.02",
            reads="4",
            assembly_support_bases="0",
        ),
    ]
    config = dict(
        expected_libraries=3, target_id="BAC00019", previously_flagged_sample="a"
    )
    return metadata, bacteria, human, config


def test_join_keeps_accepted_evidence_and_actual_prefix_denominators():
    rows = ANALYSIS.join_samples(*fixtures())
    assert rows[0]["human_compatible_pairs_per_million_raw_prefix_pairs"] == 5000
    assert rows[0]["human_compatible_pairs_per_million_qc_prefix_pairs"] == 10000
    assert rows[1]["human_compatible_pairs"] == 0
    assert rows[1]["cutibacterium_accepted"] == "false"
    assert rows[1]["cutibacterium_reads_if_accepted"] == ""
    assert rows[0]["cutibacterium_reads_if_accepted"] == "40"


def test_equal_size_wrong_cohort_and_duplicate_keys_are_rejected():
    metadata, bacteria, human, config = fixtures()
    human[2]["sample_id"] = "wrong"
    with pytest.raises(ValueError, match="full metadata cohort"):
        ANALYSIS.join_samples(metadata, bacteria, human, config)
    human[2]["sample_id"] = "a"
    with pytest.raises(ValueError, match="Nonunique"):
        ANALYSIS.join_samples(metadata, bacteria, human, config)


def test_impossible_counts_and_nonprefix_audits_are_rejected():
    metadata, bacteria, human, config = fixtures()
    human[0]["distinct_sequence_pairs"] = "11"
    with pytest.raises(ValueError, match="Invalid selected-prefix"):
        ANALYSIS.join_samples(metadata, bacteria, human, config)
    human[0]["distinct_sequence_pairs"] = "10"
    human[0]["selection"] = "random"
    with pytest.raises(ValueError, match="bounded per-source-prefix"):
        ANALYSIS.join_samples(metadata, bacteria, human, config)


def test_flagged_omission_preserves_nondetections_and_denominators():
    *_, config = fixtures()
    rows = ANALYSIS.join_samples(*fixtures())
    summary = ANALYSIS.summarize(rows, config)
    pooled = [
        r
        for r in summary
        if r["scope"] == "omit_previously_flagged" and r["study"] == "all_studies"
    ]
    assert len(pooled) == 1
    assert pooled[0]["libraries"] == 2
    assert pooled[0]["libraries_without_human_compatible_pairs"] == 1
    assert pooled[0]["raw_prefix_pairs_total"] == 4000
    assert pooled[0]["human_compatible_pairs_median"] == 0.5
