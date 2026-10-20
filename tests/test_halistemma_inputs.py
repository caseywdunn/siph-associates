"""Consequential checks for the new host reference and two-library mapping gate."""

import copy
import gzip
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name):
    spec = importlib.util.spec_from_file_location(
        name, ROOT / f"scripts/halistemma/{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


review = load("review_inputs")
reference = load("fetch_reference")


def test_reference_rejects_duplicate_ids_and_invalid_sequence(tmp_path):
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr1\nACGTNN\n>mt\nacgtry\n")
    assert reference.sequence_lengths(fasta) == {"chr1": 6, "mt": 6}
    for content in (">chr1\nACGT\n>chr1\nACGT\n", ">chr1\nACGTX\n", ">chr1\n"):
        fasta.write_text(content)
        with pytest.raises(ValueError):
            reference.sequence_lengths(fasta)


def test_pilot_mates_and_quality_must_agree(tmp_path):
    paths = [tmp_path / f"R{mate}.fq.gz" for mate in (1, 2)]
    for mate, path in enumerate(paths, 1):
        with gzip.open(path, "wt") as handle:
            handle.write(f"@pair/{mate}\nACGT\n+\nIIII\n")
    assert review.paired_count(*paths) == 1
    for content in ("@different/2\nACGT\n+\nIIII\n", "@pair/2\nACGT\n+\nIII\n"):
        with gzip.open(paths[1], "wt") as handle:
            handle.write(content)
        with pytest.raises(ValueError):
            review.paired_count(*paths)


def test_trim_rejects_changed_inputs_and_lost_pairs():
    old = {
        "sample_id": "fixture",
        "raw_pairs_expected": 10,
        "pairs_entering_fastp": 10,
        "pairs_after_fastp": 9,
        "inputs": [{"path": "R1", "size_bytes": 100, "sha256": "abc"}],
    }
    assert review.validate_trim(old, old)["status"] == "PASS"
    changed = copy.deepcopy(old)
    changed["inputs"][0]["sha256"] = "def"
    with pytest.raises(ValueError, match="checksums"):
        review.validate_trim(changed, old)
    changed = dict(old, pairs_after_fastp=8)
    with pytest.raises(ValueError, match="pairs_after_fastp"):
        review.validate_trim(changed, old)


def test_mapping_counts_and_gate():
    # Three mapped reads in three pairs: one pair maps both mates, one only one.
    # Pair-based host removal must therefore be 2/3, not the read-based 1/2.
    first = review.pilot_metrics("first", 3, 6, 3, 2)
    assert first["host_removed_pair_fraction"] == pytest.approx(2 / 3)
    assert first["mapped_read_fraction"] == 0.5
    second = review.pilot_metrics("second", 3, 6, 4, 2)
    assert (
        review.validate_pilots([first, second], {"first", "second"}, 0.5)["status"]
        == "PASS"
    )
    with pytest.raises(ValueError, match="exactly"):
        review.validate_pilots([first, first], {"first", "second"}, 0.5)
    poor = review.pilot_metrics("second", 3, 6, 0, 6)
    with pytest.raises(ValueError, match="Low host mapping"):
        review.validate_pilots([first, poor], {"first", "second"}, 0.5)
    with pytest.raises(ValueError, match="disagree"):
        review.pilot_metrics("first", 3, 5, 3, 2)
