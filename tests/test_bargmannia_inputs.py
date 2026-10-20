"""Checks for linking the new study to the public raw files and read accounting."""

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_public_upload_checksum_gate(tmp_path):
    rows = []
    for mate in (1, 2):
        path = tmp_path / f"R{mate}.fastq.gz"
        content = f"fixture mate {mate}".encode()
        path.write_bytes(content)
        rows.append(
            {
                "analysis_sample_id": "fixture",
                "mate": mate,
                "path": str(path),
                "bytes": len(content),
                "ncbi_original_md5": hashlib.md5(content).hexdigest(),
            }
        )
    manifest = tmp_path / "files.tsv"
    with manifest.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    command = [
        sys.executable,
        str(ROOT / "scripts/bargmannia/verify_raw.py"),
        "--files",
        str(manifest),
        "--sample",
        "fixture",
        "--output",
        str(tmp_path / "verified.json"),
    ]
    subprocess.run(command, check=True, capture_output=True)
    assert json.loads((tmp_path / "verified.json").read_text())["status"] == "PASS"
    # Same byte length but different content must fail before read processing.
    Path(rows[0]["path"]).write_bytes(b"x" * int(rows[0]["bytes"]))
    assert subprocess.run(command, check=False, capture_output=True).returncode != 0


def test_published_species_and_exact_selection():
    with (ROOT / "config/bargmannia_samples.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_event = {r["collection_id"]: r for r in rows}
    assert len(rows) == len(by_event) == 8
    assert by_event["D725-D8"]["species_current"] == "Bargmannia marlena"
    assert {by_event[x]["species_current"] for x in ["D611-D6", "T411-SS2"]} == {
        "Bargmannia pyrosa"
    }
    assert all(r["tissue"] for r in rows)
    assert {r["sra_run"] for r in rows} == {
        f"SRR{i}" for i in range(36270048, 36270056)
    }
    assert by_event["T676-SS3"]["preservation"] == "formalin"
    assert by_event["T676-SS3"]["voucher_id"] == ""


def test_fish_inventory_matches_consumer_columns(tmp_path):
    """The shared fish tools consume r1_path/r2_path, not r1/r2."""
    sample = "Manko2026__fixture"
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "analysis_sample_id\tspecies_current\tread_pairs\n"
        f"{sample}\tBargmannia pyrosa\t1\n"
    )
    provenance = tmp_path / "provenance/trim" / f"{sample}.json"
    provenance.parent.mkdir(parents=True)
    provenance.write_text(
        json.dumps({"parameters": {"cap_pairs": 200000000}, "pairs_after_fastp": 1})
    )
    report = tmp_path / "trim" / f"{sample}.fastp.json"
    report.parent.mkdir()
    report.write_text(
        json.dumps(
            {
                "summary": {
                    key: {"total_reads": 2}
                    for key in ["before_filtering", "after_filtering"]
                }
            }
        )
    )
    reads = tmp_path / "retained/trimmed"
    reads.mkdir(parents=True)
    for mate in (1, 2):
        (reads / f"{sample}_R{mate}.fastq.gz").write_bytes(b"fixture")
    output = tmp_path / "inventory.tsv"
    command = [
        sys.executable,
        str(ROOT / "scripts/bargmannia/read_inventory.py"),
        "--samples",
        str(samples),
        "--sample",
        sample,
        "--root",
        str(tmp_path),
        "--output",
        str(output),
    ]
    subprocess.run(command, check=True, capture_output=True)
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1 and rows[0]["read_pair_id"] == "capped_trimmed"
    # These are the exact columns used by both independent fish-screen consumers.
    for mate in (1, 2):
        path = Path(rows[0][f"r{mate}_path"])
        assert path.stat().st_size == int(rows[0][f"r{mate}_bytes"])
    assert int(rows[0]["pairs"]) == 1
