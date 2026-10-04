#!/usr/bin/env python3
"""Regressions for overlapping eukaryote evidence and honest taxonomic naming.

Run directly with Python, as with the repository's existing regression scripts.
Only small temporary TSV fixtures are used; no reference databases are needed.
"""

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
from eukaryote_reporting import annotate_detections  # noqa: E402


def record(sample: str, unit: str, grade: str = "validated") -> dict:
    return {"sample_id": sample, "reporting_unit": unit, "grade": grade}


def test_overlapping_units() -> None:
    source = [
        record("one", "algae"),
        record("one", "algae;red"),
        record("one", "algae;red;group_a", "high_confidence"),
        record("one", "algae;green"),
        record("one", "algae_other"),
        record("two", "algae"),
        record("two", "algae;red", "trace"),
    ]
    result = annotate_detections(source)
    assert len(result) == len(source), "ambiguous evidence must be retained"
    assert [r["count_as_detection"] for r in result] == [
        "false", "false", "true", "true", "true", "true", "false"
    ]
    assert result[0]["overlapping_descendant_units"] == "algae;green|algae;red|algae;red;group_a"
    assert result[0]["reporting_status"] == "unresolved_ancestor"
    assert result[5]["reporting_status"] == "detection", "trace descendants cannot suppress presence"
    assert "count_as_detection" not in source[0], "source evidence must not be mutated"
    assert annotate_detections([]) == []


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_grade_and_validation() -> None:
    with tempfile.TemporaryDirectory(prefix="eukaryote-reporting-") as temporary:
        root = Path(temporary)
        config = root / "config.json"
        settings = {
            "presence": {
                "validated_reads": {"min_identity_band": 97, "min_read_pairs": 37},
                "trace_reads": {"min_identity_band": 97, "min_read_pairs": 2},
            },
            "role_field": {"parasite": [], "prey": ["Copepoda"]},
        }
        config.write_text(json.dumps(settings))
        parent = "Eukaryota;Archaeplastida;Rhodophyceae;Florideophycidae"
        child = parent + ";Rhodymeniophycidae"
        order = "Eukaryota;Amorphea;Metazoa;Arthropoda;Crustacea;Copepoda;Calanoida"
        reads = root / "reads.tsv"
        write_tsv(
            reads,
            ["sample_id", "lineage", "min_identity_band", "read_pairs"],
            [
                {"sample_id": sample, "lineage": lineage, "min_identity_band": 99, "read_pairs": pairs}
                for sample, lineage, pairs in [
                    ("one", parent, 324), ("one", child, 961),
                    ("two", parent, 40), ("two", child, 5), ("three", order, 45),
                ]
            ],
        )
        assembled = root / "assembled.tsv"
        write_tsv(
            assembled,
            ["sample_id", "source", "sequence_id", "lineage", "identity", "aligned_bases"],
            [{"sample_id": "one", "source": "phyloflash", "sequence_id": "seq1",
              "lineage": child, "identity": 98, "aligned_bases": 1500}],
        )
        verification = root / "verification.tsv"
        write_tsv(
            verification,
            ["sequence_id", "rank", "organism", "title", "subject", "identity"],
            [{"sequence_id": "seq1", "rank": 1, "organism": "red alga", "title": "SSU",
              "subject": "test-accession", "identity": 98}],
        )
        grades = root / "grades.tsv"
        common = ["--config", str(config), "--assembled", str(assembled), "--verification", str(verification)]
        subprocess.run(
            [sys.executable, "-B", str(ROOT / "workflow/scripts/grade_eukaryotes.py"),
             *common, "--reads", str(reads), "--output", str(grades)],
            check=True, capture_output=True, text=True,
        )
        with grades.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 5
        assert sum(r["count_as_detection"] == "true" for r in rows) == 3
        ancestor = next(r for r in rows if r["sample_id"] == "one" and r["reporting_unit"] == parent[10:])
        assert ancestor["read_pairs_ge97"] == "324", "do not allocate ambiguous reads to a descendant"
        assert ancestor["grade"] == "validated" and ancestor["count_as_detection"] == "false"
        descendant = next(r for r in rows if r["grade"] == "high_confidence")
        assert descendant["read_pairs_ge97"] == "961"
        copepod = next(r for r in rows if r["sample_id"] == "three")
        assert copepod["lca_terminal_taxon"] == "Calanoida"
        assert copepod["named_lineage"] == order
        assert copepod["naming_ceiling"] == "genus" and "naming_depth" not in copepod
        done = root / "done.tsv"
        validate = [
            sys.executable, "-B", str(ROOT / "workflow/scripts/validate_eukaryote_gate.py"),
            *common, "--grades", str(grades), "--output", str(done),
        ]
        subprocess.run(validate, check=True, capture_output=True, text=True)
        summary = dict(line.split("\t") for line in done.read_text().splitlines())
        assert summary["validated_evidence_rows"] == "4"
        assert summary["validated_detections"] == "3"
        assert summary["validated_ancestor_evidence_rows"] == "1"
        ancestor["count_as_detection"] = "true"
        write_tsv(grades, list(rows[0]), rows)
        rejected = subprocess.run(validate, check=False, capture_output=True, text=True)
        assert rejected.returncode != 0
        assert "incorrect detection count flag" in rejected.stdout


if __name__ == "__main__":
    test_overlapping_units()
    test_grade_and_validation()
    print("status\tPASS")
