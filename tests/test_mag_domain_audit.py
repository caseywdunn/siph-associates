"""Guard decisions that could turn mixed or unmapped evidence into an assignment."""

import importlib.util
import sys
from pathlib import Path

folder = Path(__file__).resolve().parents[1] / "scripts/mag_domain_audit"
sys.path.insert(0, str(folder))
spec = importlib.util.spec_from_file_location("audit_classify", folder / "classify.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def hit(label: str, score: float = 100, length: int = 100) -> dict:
    return {"domains": {label}, "bitscore": score, "length": length, "evalue": 1e-20}


def classify(hits: list[dict]) -> str:
    return audit.classify_hits(hits, 0.95, 50, 1e-5)[0]


def test_cross_domain_near_tie_stays_ambiguous():
    assert classify([hit("Archaea"), hit("Eukaryota", 96)]) == "ambiguous"


def test_stronger_eukaryote_match_wins():
    assert classify([hit("Eukaryota"), hit("Archaea", 80)]) == "Eukaryota"


def test_unknown_taxonomy_is_not_ignored():
    assert (
        classify([hit("Archaea"), hit("unresolved_taxonomy", 98)])
        == "unresolved_taxonomy"
    )


def test_short_hits_do_not_establish_domain():
    assert classify([hit("Archaea", length=35)]) == "no_qualifying_hit"
    assert classify([]) == "no_qualifying_hit"


def test_taxonomy_walk_includes_eukaryotes_and_unknown_ids():
    parents = {123: 2759, 2759: 131567, 131567: 1, 1: 1}
    assert audit.domain(123, parents) == "Eukaryota"
    assert audit.domain(999, parents) == "unresolved_taxonomy"
