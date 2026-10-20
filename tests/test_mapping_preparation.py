"""QC preparation must preserve mapping holds and reject ambiguous source records."""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
spec = importlib.util.spec_from_file_location(
    "preparation", ROOT / "scripts/prepare_mapping_review.py"
)
preparation = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preparation)


def test_mapping_hold_cannot_be_implicitly_lifted(tmp_path):
    path = tmp_path / "hold.json"
    record = {
        "status": "hold",
        "automatic_remapping": False,
        "requires_explicit_mapping_authorization": True,
    }
    path.write_text(json.dumps(record))
    preparation.check_hold(str(path))
    for key, value in [
        ("status", "ready"),
        ("automatic_remapping", True),
        ("requires_explicit_mapping_authorization", False),
    ]:
        path.write_text(json.dumps(dict(record, **{key: value})))
        with pytest.raises(ValueError, match="mapping hold"):
            preparation.check_hold(str(path))


def test_candidate_contig_sequences_and_ids_are_validated(tmp_path):
    path = tmp_path / "bin.fa"
    path.write_text(">sample__contig1\nACGTNN\n>sample__contig2\nGGCC\n")
    rows = preparation.fasta_metrics(path)
    assert rows[0]["length"] == 6 and rows[0]["gc_fraction"] == 0.5
    assert rows[0]["n_fraction"] == pytest.approx(1 / 3)
    assert rows[1]["gc_fraction"] == 1
    for data in ("", ">x\nACGT\n>x\nACGT\n", ">x\nACGTZ\n", "ACGT\n"):
        path.write_text(data)
        with pytest.raises(ValueError):
            preparation.fasta_metrics(path)


def test_full_mapping_never_substitutes_assembly_subsample(tmp_path):
    sample = "fixture__not_real"
    provenance = tmp_path / "provenance/trim" / f"{sample}.json"
    provenance.parent.mkdir(parents=True)
    provenance.write_text(
        json.dumps({"sample_id": sample, "pairs_after_fastp": 190000000})
    )
    # A retained 25M-pair assembly input is present but is not a full mapping input.
    assembly = tmp_path / "phase3_cohort/inputs"
    assembly.mkdir(parents=True)
    for mate in (1, 2):
        (assembly / f"{sample}_R{mate}.fastq.gz").write_bytes(b"fixture")
    row = {
        "sample_id": sample,
        "host_route": "none",
        "input_root": str(tmp_path),
        "retained_pairs": 25000000,
        "species_current": "Host",
        "r1_paths": "raw1",
        "r2_paths": "raw2",
        "raw_pairs": 220000000,
    }
    record = preparation.mapping_input(row, {})
    assert record["read_pairs"] == 190000000
    assert (
        record["availability"] == "regenerate_trimmed_reads_when_mapping_is_authorized"
    )
    assert record["r1"] == record["r2"] == ""
