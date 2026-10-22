"""Known-answer checks for the stringent human-read decision."""

import importlib.util
from pathlib import Path

import pysam
import pytest

SPEC = importlib.util.spec_from_file_location(
    "human_filter",
    Path(__file__).parents[1] / "scripts/human_screen/filter_alignments.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def read(flag, nm=1, mapq=42, cigar="150M"):
    record = pysam.AlignedSegment()
    record.query_name = "pair"
    record.query_sequence = "ACGT" * 37 + "AC"
    record.flag = flag
    record.reference_id = 0
    record.reference_start = 100 if flag == 99 else 300
    record.mapping_quality = mapq
    record.cigarstring = cigar
    record.next_reference_id = 0
    record.template_length = 350 if flag == 99 else -350
    record.set_tag("NM", nm)
    return record


def test_near_exact_proper_pair_passes():
    assert MODULE.passes_pair([read(99), read(147)])


@pytest.mark.parametrize(
    "change",
    [
        "divergent",
        "low_mapq",
        "clipped",
        "unpaired",
        "short",
        "low_complexity",
        "one_mate",
    ],
)
def test_nonhuman_like_or_ambiguous_pairs_fail(change):
    pair = [read(99), read(147)]
    if change == "divergent":
        pair[1].set_tag("NM", 2)
    elif change == "low_mapq":
        pair[1].mapping_quality = 20
    elif change == "clipped":
        pair[1].cigarstring = "10S140M"
    elif change == "unpaired":
        pair[1].flag = 145
    elif change == "short":
        pair[1].query_sequence = "ACGT" * 20
        pair[1].cigarstring = "80M"
    elif change == "low_complexity":
        pair[1].query_sequence = "A" * 150
    elif change == "one_mate":
        pair = pair[:1]
    assert not MODULE.passes_pair(pair)


NT_SPEC = importlib.util.spec_from_file_location(
    "human_nt", Path(__file__).parents[1] / "scripts/human_screen/review_nt.py"
)
NT = importlib.util.module_from_spec(NT_SPEC)
NT_SPEC.loader.exec_module(NT)


def hit(taxon, score=200):
    return {
        "pident": "100",
        "qstart": "1",
        "qend": "150",
        "qlen": "150",
        "bitscore": str(score),
        "staxids": taxon,
        "sscinames": "taxon",
    }


def test_truncated_human_only_hits_remain_unresolved():
    result = NT.classify_query([hit("9606")], 100)
    assert result["query_classification"] == "unresolved_truncation_or_unknown_taxonomy"


def test_unknown_taxonomy_cannot_hide_alternative():
    result = NT.classify_query([hit("9606"), hit("N/A")], 2)
    assert result["query_classification"] == "unresolved_truncation_or_unknown_taxonomy"


def test_nonhuman_near_tie_is_preserved():
    result = NT.classify_query([hit("9606"), hit("9606", 199), hit("7955", 195)], 3)
    assert result["query_classification"] == "human_with_alternatives"


def test_human_pair_requires_both_mates():
    human = NT.classify_query([hit("9606")], 1)
    nonhuman = NT.classify_query([hit("7955")], 1)
    assert NT.classify_pair([human, human]) == "both_mates_human_only_near_best"
    assert NT.classify_pair([human, nonhuman]) == "conflicting_mates"
    assert (
        NT.classify_pair([human, NT.classify_query([], 0)]) == "incomplete_confirmation"
    )
