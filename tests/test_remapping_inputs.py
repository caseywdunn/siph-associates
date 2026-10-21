"""Mapping checks must reject changed read sources despite matching filenames."""

import copy
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from verify_remapping_inputs import (
    verify_authorization,
    verify_read_record,
    verify_regeneration,
)


@pytest.fixture
def reads(tmp_path):
    outputs = []
    row = {
        "sample_id": "fixture",
        "host_route": "none",
        "read_pairs": "100",
        "cap_pairs": "200000000",
        "cap_seed": "20260829",
    }
    for mate in (1, 2):
        path = tmp_path / f"fixture_R{mate}.fastq.gz"
        path.write_bytes(b"fixture")
        row[f"r{mate}"] = str(path)
        outputs.append({"path": str(path), "size_bytes": 7})
    record = {
        "sample_id": "fixture",
        "software": "fastp 0.23.2",
        "outputs": outputs,
        "inputs": [{"path": "raw.fastq.gz", "size_bytes": 90, "sha256": "original"}],
        "raw_pairs_expected": 110,
        "pairs_entering_fastp": 110,
        "pairs_after_fastp": 100,
        "parameters": {"cap_pairs": 200000000, "cap_seed": 20260829},
    }
    return row, record


def test_regeneration_rejects_changed_raw_sequence_with_same_counts(reads):
    row, record = reads
    current = copy.deepcopy(record)
    verify_regeneration(row, record, current)
    current["inputs"][0]["sha256"] = "different"
    with pytest.raises(ValueError, match="Raw input identities changed"):
        verify_regeneration(row, record, current)


def test_full_read_check_rejects_assembly_subset(reads):
    row, record = reads
    record["pairs_after_fastp"] = 25
    with pytest.raises(ValueError, match="Full trimmed pair count mismatch"):
        verify_read_record(row, record)


def test_host_read_check_rejects_full_read_source(reads):
    row, record = reads
    row["host_route"] = "H_rubrum"
    with pytest.raises(ValueError, match="Host-depleted read source mismatch"):
        verify_read_record(row, record)


def test_regeneration_rejects_different_seed(reads):
    row, record = reads
    current = copy.deepcopy(record)
    current["parameters"]["cap_seed"] += 1
    with pytest.raises(ValueError, match="capping differs"):
        verify_regeneration(row, record, current)


def test_preflight_rejects_hold_without_reading_catalog(tmp_path):
    authorization = tmp_path / "authorization.json"
    authorization.write_text(
        json.dumps({"status": "hold", "automatic_remapping": False})
    )
    with pytest.raises(ValueError, match="explicit user authorization"):
        verify_authorization({"authorization": str(authorization)})
