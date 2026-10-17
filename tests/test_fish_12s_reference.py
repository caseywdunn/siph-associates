"""Known-answer tests for strand orientation and fixed-base substitutions."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/fish_12s"))
from assess_reference import MOTIF, closest_window, reverse_complement


def test_reverse_complement_matches_pilot() -> None:
    assert reverse_complement(MOTIF) == "GCATA..GGGGTAT..AATCC..GTTTGT"


def test_wildcards_are_not_substitutions() -> None:
    seq = MOTIF.replace(".", "A")
    assert closest_window(seq)[0] == 0
    assert closest_window(reverse_complement(seq))[0] == 0
    assert closest_window(reverse_complement(seq))[2] == "-"


def test_cartilaginous_variant_has_two_fixed_base_changes() -> None:
    seq = "CCAAACTAGGATTAGATACCCTACTATGC"
    assert closest_window(seq)[0] == 2


def test_flanking_sequence_coordinates() -> None:
    seq = "TTTTTTT" + MOTIF.replace(".", "A") + "TTTT"
    distance, position, strand, window = closest_window(seq)
    assert (distance, position, strand) == (0, 7, "+")
    assert len(window) == 29
