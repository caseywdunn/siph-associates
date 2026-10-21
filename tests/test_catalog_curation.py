"""Protect viral coordinates and historical catalog identifiers during curation."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from catalog_ids import stable_ids
from prepare_catalog_curation import viral_regions


def test_retired_ids_are_not_reused():
    previous = {"MAGSP0001": {"a"}, "MAGSP0002": {"removed"}, "MAGSP0003": {"b"}}
    assert stable_ids([{"b"}, {"new"}, {"a"}], previous, "MAGSP", 4) == [
        "MAGSP0003",
        "MAGSP0004",
        "MAGSP0001",
    ]


def test_largest_split_keeps_old_id_and_merge_gets_one_unique_id():
    previous = {"VOTU00001": {"a", "b", "c"}, "VOTU00002": {"d"}}
    assert stable_ids([{"a"}, {"b", "c", "d"}], previous, "VOTU", 5) == [
        "VOTU00003",
        "VOTU00001",
    ]


def test_separated_viral_regions_are_not_concatenated():
    row = {
        "contig_id": "test",
        "contig_length": "100",
        "provirus": "Yes",
        "region_types": "viral,host,viral",
        "region_coords_bp": "1-30,31-70,71-100",
    }
    assert viral_regions(row) == [(1, 30), (71, 100)]
    row["region_coords_bp"] = "1-30,32-70,71-100"
    with pytest.raises(ValueError, match="partition"):
        viral_regions(row)


def test_untrimmed_sequence_keeps_complete_inclusive_interval():
    assert viral_regions({"contig_length": "100", "provirus": "No"}) == [(1, 100)]
