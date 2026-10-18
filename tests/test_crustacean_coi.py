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


def test_long_hit_list_only_ambiguous_when_near_best_is_truncated():
    metadata = {str(i): {"target": "1", "group": "Copepoda"} for i in range(200)}
    hits = [
        {"reference_id": str(i), "length": 150, "identity": 98, "bitscore": 280}
        for i in range(200)
    ]
    assert summary.classify(hits, metadata)[0] == "ambiguous"
    hits[-1]["bitscore"] = 150
    assert summary.classify(hits, metadata)[0] == "crustacean"
