"""Biological and classification checks for the isolated COI pilot."""

import importlib.util
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature


def module(name):
    path = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "crustacean_coi"
        / f"{name}.py"
    )
    spec = importlib.util.spec_from_file_location(name, path)
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


reference = module("prepare_reference")
summary = module("summarize")
verification = module("verify")


@pytest.mark.parametrize(
    "name,expected",
    [
        ("cytochrome c oxidase subunit I", True),
        ("cytochrome c oxidase subunit 1", True),
        ("cytochrome c oxidase subunit II", False),
        ("cytochrome c oxidase subunit III", False),
        ("cytochrome b", False),
    ],
)
def test_coi_annotation_excludes_other_mitochondrial_genes(name, expected):
    assert (
        reference.is_coi(SeqFeature(type="CDS", qualifiers={"product": [name]}))
        == expected
    )


def test_invertebrate_mitochondrial_code_and_reverse_frame():
    # ATA=M, TGA=W and AGA=S differ from the standard nuclear code.
    sequence = "ATATGAAGA"
    assert reference.translate_frame(sequence, 1) == "MWS"
    assert reference.translate_frame("A" + sequence, 2) == "MWS"
    assert (
        reference.translate_frame(str(Seq("A" + sequence).reverse_complement()), -2)
        == "MWS"
    )


def test_host_tied_with_crustacean_is_ambiguous():
    metadata = {
        "prey": {"target": "1", "group": "Copepoda"},
        "host": {"target": "0", "group": "Cnidaria"},
    }
    hits = [
        {"reference_id": r, "length": 150, "identity": 98, "bitscore": 280}
        for r in metadata
    ]
    assert summary.classify(hits, metadata)[0] == "ambiguous"
    hits[1]["bitscore"] = 180
    assert summary.classify(hits, metadata)[0] == "crustacean"


def test_discordant_mates_are_ambiguous():
    assert summary.combine(["crustacean", "non_target"]) == "ambiguous"
    assert summary.combine(["crustacean", "unresolved"]) == "crustacean"


def test_crustacean_lineage_includes_barnacles_but_excludes_insects():
    assert verification.origin(["Altocrustacea", "Cirripedia"]) == "crustacean"
    assert verification.origin(["Oligostraca", "Ostracoda"]) == "crustacean"
    assert verification.origin(["Altocrustacea", "Hexapoda"]) == "other_eukaryote"


def test_divergent_non_target_protein_does_not_reject_potential_prey():
    assert (
        verification.decide_pair(["other_eukaryote", "weak"], ["weak", "unresolved"])
        == "unresolved_candidate"
    )


def test_bacterial_protein_cannot_establish_crustacean_detection():
    assert (
        verification.decide_pair(["prokaryote", "prokaryote"], ["ambiguous", "weak"])
        == "bacterial_or_archaeal_oxidase"
    )
    assert (
        verification.decide_pair(["prokaryote", "crustacean"], ["weak", "crustacean"])
        == "ambiguous"
    )


def test_support_must_come_from_same_mate():
    assert (
        verification.decide_pair(["crustacean", "unresolved"], ["weak", "crustacean"])
        == "unresolved_candidate"
    )
    assert (
        verification.decide_pair(["crustacean", "unresolved"], ["crustacean", "weak"])
        == "crustacean_supported"
    )


def test_quarantined_nucleotide_reference_cannot_support_crustacean_call():
    metadata = {"suspect": {"origin": "unresolved"}, "prey": {"origin": "crustacean"}}
    alignments = [
        {"reference_id": "suspect", "identity": 100, "length": 151, "bitscore": 300},
        {"reference_id": "prey", "identity": 98, "length": 151, "bitscore": 290},
    ]
    assert (
        verification.nucleotide_evidence(
            alignments, metadata, {"read_near_best_fraction": 0.95}
        )
        == "ambiguous"
    )


def test_short_match_cannot_replace_stronger_divergent_competitor():
    metadata = {
        "strong": {"target": "1", "group": "Decapoda"},
        "weak": {"target": "1", "group": "Amphipoda"},
    }
    hits = [
        {"reference_id": "strong", "length": 150, "identity": 75.333, "bitscore": 104},
        {"reference_id": "weak", "length": 98, "identity": 80.612, "bitscore": 92.4},
    ]
    assert summary.classify(hits, metadata) == ("unresolved", ["Decapoda"], ["strong"])


def test_non_target_competitor_below_identity_cutoff_is_not_discarded():
    metadata = {
        "host": {"target": "0", "group": "Cnidaria"},
        "prey": {"target": "1", "group": "Copepoda"},
    }
    hits = [
        {"reference_id": "host", "length": 151, "identity": 79, "bitscore": 120},
        {"reference_id": "prey", "length": 100, "identity": 81, "bitscore": 119},
    ]
    assert summary.classify(hits, metadata)[0] == "ambiguous"


def test_long_hit_list_only_ambiguous_when_near_best_is_truncated():
    metadata = {str(i): {"target": "1", "group": "Copepoda"} for i in range(200)}
    hits = [
        {"reference_id": str(i), "length": 150, "identity": 98, "bitscore": 280}
        for i in range(200)
    ]
    assert summary.classify(hits, metadata)[0] == "ambiguous"
    assert summary.classify(hits, metadata, max_hits=None)[0] == "crustacean"
    hits[-1]["bitscore"] = 150
    assert summary.classify(hits, metadata)[0] == "crustacean"
