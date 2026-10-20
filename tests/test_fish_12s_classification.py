"""Regression checks for the fish/non-fish decision and final control gate."""

import csv
import importlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts/fish_12s"))
classification = importlib.import_module("classify_candidates")
metadata = importlib.import_module("reference_metadata")


def hit(fish, score=200, length=151, identity=100):
    return {
        "fish": fish,
        "bitscore": score,
        "length": length,
        "identity": identity,
        "name": str(fish),
    }


def test_competing_frog_cannot_be_called_fish():
    assert (
        classification.evidence([hit(True), hit(False)])["category"]
        == "other_fish_similarity"
    )
    fish = classification.evidence([hit(True)])
    frog = classification.evidence([hit(False)])
    assert classification.pair_category({1: fish, 2: frog}) == "other_fish_similarity"


def test_missing_mate_hit_is_allowed_but_truncated_evidence_is_not_close():
    fish = classification.evidence([hit(True)])
    assert (
        classification.pair_category({1: fish, 2: classification.evidence([])})
        == "long_close_fish_match"
    )
    assert (
        classification.evidence([hit(True)] * 200)["category"]
        == "other_fish_similarity"
    )


def test_accession_lineages_do_not_confuse_tetrapods_with_fish():
    assert metadata.is_fish(["Vertebrata", "Actinopteri", "Teleostei"])
    assert metadata.is_fish(["Vertebrata", "Chondrichthyes"])
    assert metadata.is_fish(["Vertebrata", "Sarcopterygii", "Dipnoi"])
    assert not metadata.is_fish(
        ["Vertebrata", "Sarcopterygii", "Tetrapoda", "Amphibia"]
    )


def test_missing_alignment_file_fails(tmp_path):
    with pytest.raises(FileNotFoundError):
        classification.read_hits(tmp_path / "absent.tsv", {})


def test_unknown_reference_fails(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text("query\tmissing_reference\t100\t151\n")
    with pytest.raises(KeyError):
        classification.read_hits(path, {})


def test_refseq_identifier_and_reverse_marker_overlap(tmp_path):
    path = tmp_path / "hits.tsv"
    path.write_text(
        "q\tref|NC_001573.1|\t100\t101\t0\t0\t1\t101\t200\t100\t0\t200\t151\t1000\tFrog\n"
    )
    refs = {
        "NC_001573.1": {
            "name": "Frog",
            "fish": "0",
            "marker_intervals_1based": "120-150",
        }
    }
    result = classification.read_hits(path, refs)["q"][0]
    assert result["marker_overlap_nt"] == 31
    assert not result["fish"]


@pytest.mark.parametrize("mode", ["valid", "missing", "empty", "frog_as_fish"])
def test_control_gate_checks_alignments(tmp_path, monkeypatch, mode):
    gate = importlib.import_module("check_cohort_controls")
    out = tmp_path / "union/synthetic/host_depleted"
    (out / "alignment").mkdir(parents=True)
    names = [
        "motif_R1",
        "motif_reverse_R2",
        "fish_away_from_motif",
        "cartilaginous_motif_variant",
        "frog",
    ]
    truth = tmp_path / "truth.json"
    truth.write_text(json.dumps({"expected_pair_ids": names}))
    with (out / "candidate_pairs.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["candidate_id", "read_pair_id"])
        writer.writerows((f"c{i}", name) for i, name in enumerate(names))
    (out / "candidates.fasta").write_text(
        "".join(f">c{i}_{mate}\nACGT\n" for i in range(5) for mate in [1, 2])
    )
    (out / "run.json").write_text(
        json.dumps(
            {
                "input_pairs": 9,
                "method_pairs": {
                    "regex": 3,
                    "bbduk25": 5,
                    "sortmerna": 5,
                    "bowtie2": 5,
                },
            }
        )
    )
    refs = tmp_path / "refs.tsv"
    refs.write_text(
        "reference_id\tname\tfish\tmarker_intervals_1based\nfish\tFish\t1\t1-1000\nfrog\tFrog\t0\t1-1000\n"
    )
    if mode != "missing":
        rows = []
        for i in range(5):
            ref = "frog" if i == 4 and mode != "frog_as_fish" else "fish"
            rows.append(
                f"c{i}_1\t{ref}\t100\t151\t0\t0\t1\t151\t1\t151\t0\t200\t151\t1000\tTitle\n"
            )
        (out / "alignment/hits.tsv").write_text(
            "" if mode == "empty" else "".join(rows)
        )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "check",
            "--root",
            str(tmp_path),
            "--truth",
            str(truth),
            "--reference-metadata",
            str(refs),
        ],
    )
    if mode == "valid":
        gate.main()
        assert json.loads((tmp_path / "CONTROLS_PASSED").read_text())[
            "competitive_classification_checked"
        ]
    else:
        with pytest.raises((FileNotFoundError, ValueError)):
            gate.main()
        assert not (tmp_path / "CONTROLS_PASSED").exists()
