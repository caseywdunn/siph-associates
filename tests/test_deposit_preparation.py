"""Consequential checks for reviewed MAG edits and original coordinate handling."""

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/prepare_deposit_mag.py"
SPEC = importlib.util.spec_from_file_location("prepare_deposit_mag", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_terminal_trims_use_original_one_based_inclusive_coordinates():
    sequence = "A" * 500
    trims = [
        {"contig_id": "left", "start_1based": "1", "end_1based": "27",
         "original_length": "500", "reason": "adaptor"},
        {"contig_id": "right", "start_1based": "482", "end_1based": "500",
         "original_length": "500", "reason": "adaptor"},
    ]
    edited, events = MODULE.curate_mag({"sample__left": sequence, "sample__right": sequence}, [], trims, [])
    assert len(edited["sample__left"]) == 473
    assert len(edited["sample__right"]) == 481
    assert [r["retained_length"] for r in events] == [473, 481]


def test_palindrome_repair_requires_exact_duplicate_reverse_complement():
    half = "AACCGGTTAC" * 30
    palindrome = half + MODULE.reverse_complement(half)
    repair = {"contig_id": "palindrome", "action": "retain_first_half_of_exact_reverse_complement_palindrome",
              "original_length": "600", "retained_length": "300", "reason": "reviewed palindrome"}
    edited, _ = MODULE.curate_mag({"sample__palindrome": palindrome}, [], [], [repair])
    assert edited["sample__palindrome"] == half
    with pytest.raises(ValueError, match="not an exact palindrome"):
        MODULE.curate_mag({"sample__palindrome": "C" + palindrome[1:]}, [], [], [repair])


def test_removal_and_trim_cannot_both_target_same_contig():
    removals = [{"contig_id": "x", "length": "500", "reason": "misbin"}]
    trims = [{"contig_id": "x", "start_1based": "1", "end_1based": "27",
              "original_length": "500", "reason": "adaptor"}]
    with pytest.raises(ValueError, match="multiply edited"):
        MODULE.curate_mag({"x": "A" * 500, "y": "C" * 500}, removals, trims, [])


def test_internal_adaptor_requires_an_explicit_new_decision():
    trim = {"contig_id": "x", "start_1based": "100", "end_1based": "120",
            "original_length": "500", "reason": "adaptor"}
    with pytest.raises(ValueError, match="terminal adaptor"):
        MODULE.curate_mag({"x": "A" * 500}, [], [trim], [])


def test_stale_reviewed_length_is_not_silently_accepted():
    removal = {"contig_id": "x", "length": "501", "reason": "misbin"}
    with pytest.raises(ValueError, match="Reviewed length changed"):
        MODULE.curate_mag({"x": "A" * 500, "y": "C" * 500}, [removal], [], [])
