"""Checks for alignment-coordinate transformations affecting taxon comparisons."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parents[1] / "scripts" / "waki_comparison"))
from compare_waki_assemblies import subject_map
from summarize_waki_comparison import comparison_taxon, shared_comparison


def hit(qseq, sseq, start, end):
    return {"qseq": qseq, "sseq": sseq, "sstart": str(start), "send": str(end)}


def test_forward_subject_coordinates_retain_subject_insertions():
    assert subject_map(hit("AC-G", "ACTG", 4, 7)) == {
        4: ("A", "A"),
        5: ("C", "C"),
        6: ("-", "T"),
        7: ("G", "G"),
    }


def test_reverse_coordinates_and_bases_are_normalized_together():
    assert subject_map(hit("AC-G", "ACTG", 7, 4)) == {
        7: ("T", "T"),
        6: ("G", "G"),
        5: ("-", "A"),
        4: ("C", "C"),
    }


def test_reference_insertions_do_not_have_subject_coordinates():
    assert subject_map(hit("ACGT", "AC-T", 1, 3)) == {
        1: ("A", "A"),
        2: ("C", "C"),
        3: ("T", "T"),
    }


def test_inconsistent_subject_coordinates_fail():
    with pytest.raises(ValueError, match="coordinates disagree"):
        subject_map(hit("ACGT", "ACGT", 1, 3))


def test_shared_interval_removes_nonoverlapping_flanks():
    scores, positions = shared_comparison(
        [
            hit("AACCGG", "AACCGG", 1, 6),
            hit("CTGGAA", "CCGGAA", 3, 8),
        ]
    )
    assert positions == [3, 4, 5, 6]
    assert [row["identity"] for row in scores] == [100, 75]


def test_reverse_and_forward_alignments_share_identical_coordinates():
    scores, positions = shared_comparison(
        [
            hit("ACGT", "ACGT", 1, 4),
            hit("AC-T", "ACGT", 4, 1),
        ]
    )
    assert positions == [1, 2, 3, 4]
    assert [row["identity"] for row in scores] == [100, 75]


def test_reference_aliases_do_not_create_independent_competitors():
    assert comparison_taxon({"paper_taxon": "Hemiuridae sp. A"}) == "Dinurus barbatus"
    assert (
        comparison_taxon({"paper_taxon": "Prodistomum Type 4"}) == "Prodistomum Type 4"
    )
