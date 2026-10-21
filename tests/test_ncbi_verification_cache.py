"""Query identity, not the library name, controls reuse of dated NCBI evidence."""

import csv
import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "ncbi_cache",
    Path(__file__).resolve().parents[1] / "workflow/scripts/ncbi_verification_cache.py",
)
cache = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cache)


def test_identical_sequence_reuses_hits_and_retains_original_date(tmp_path):
    fasta = tmp_path / "old.fasta"
    fasta.write_text(">old\nACGT\n>nohit\nAAAA\n")
    table = tmp_path / "old.tsv"
    with table.open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sequence_id", "rank", "query_date", "rid", "subject"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(
            [
                dict(
                    sequence_id="old",
                    rank=1,
                    query_date="2026-09-30",
                    rid="dated",
                    subject="A",
                ),
                dict(
                    sequence_id="old",
                    rank=2,
                    query_date="2026-09-30",
                    rid="dated",
                    subject="B",
                ),
                dict(
                    sequence_id="nohit",
                    rank=0,
                    query_date="2026-09-30",
                    rid="empty",
                    subject="no significant hit",
                ),
            ]
        )
    reused, pending = cache.cached_verifications(
        {"new_name": "ACGT", "old": "ACGA", "other": "AAAA"}, fasta, table
    )
    assert pending == ["old"]
    assert [r["sequence_id"] for r in reused] == ["new_name", "new_name", "other"]
    assert all(r["query_date"] == "2026-09-30" for r in reused)
    assert reused[0]["rid"] == "dated"
    assert reused[0]["verification_source"].endswith("#old")
    assert reused[2]["rank"] == "0"


def test_cache_only_cli_never_opens_network(tmp_path, monkeypatch):
    import runpy
    import sys
    import urllib.request

    import pytest

    fasta = tmp_path / "old.fasta"
    fasta.write_text(">old\nACGT\n")
    table = tmp_path / "old.tsv"
    table.write_text(
        "sequence_id\trank\tquery_date\trid\tsubject\nold\t1\t2026-09-30\tdated\tA\n"
    )
    query = tmp_path / "query.fasta"
    output = tmp_path / "verification.tsv"
    script = (
        Path(__file__).resolve().parents[1]
        / "workflow/scripts/verify_eukaryotes_ncbi.py"
    )
    monkeypatch.syspath_prepend(str(script.parent))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(script),
            "--fasta",
            str(query),
            "--email",
            "test@example.invalid",
            "--output",
            str(output),
            "--cache-fasta",
            str(fasta),
            "--cache-table",
            str(table),
            "--cache-only",
        ],
    )

    def no_network(*args, **kwargs):
        raise AssertionError("Cache-only verification must not open a network request")

    monkeypatch.setattr(urllib.request, "urlopen", no_network)
    query.write_text(">new_name\nACGT\n")
    runpy.run_path(str(script), run_name="__main__")
    result = list(csv.DictReader(output.open(), delimiter="\t"))
    assert result[0]["sequence_id"] == "new_name"
    assert result[0]["query_date"] == "2026-09-30"
    output.unlink()
    query.write_text(">changed\nACGA\n")
    with pytest.raises(SystemExit, match="Cache-only verification"):
        runpy.run_path(str(script), run_name="__main__")
    assert not output.exists()
