"""A replacement input must supersede, rather than duplicate, a library screen."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/fish_12s"))
from finalize_report import select_full_read_libraries  # noqa: E402


def test_depleted_replacements_remove_only_corresponding_full_read_libraries():
    original = [{"analysis_sample_id": s} for s in ["CWD12", "CWD17", "other"]]
    replacements = [{"analysis_sample_id": s} for s in ["CWD12", "CWD17"]]
    assert select_full_read_libraries(original, replacements) == [original[2]]
    assert len(original) == 3


@pytest.mark.parametrize("ids", [["missing"], ["CWD12", "CWD12"]])
def test_replacement_identifiers_must_match_once(ids):
    with pytest.raises(ValueError):
        select_full_read_libraries(
            [{"analysis_sample_id": "CWD12"}],
            [{"analysis_sample_id": s} for s in ids],
        )
