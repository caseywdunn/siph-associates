"""Coverage reports must distinguish reused full reads from regenerated reads."""

import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def write_table(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


@pytest.mark.parametrize(
    "availability,expected",
    [("ready", "trimmed_retained"), ("missing", "trimmed_regenerated")],
)
def test_aggregate_and_validate_retained_or_regenerated_input(
    tmp_path, availability, expected
):
    sample = "fixture"
    results, provenance = tmp_path / "results", tmp_path / "provenance"
    for kind, target, heading in [
        ("bacterial", "BAC00001", "Genome"),
        ("viral", "VOTU00001", "Contig"),
    ]:
        write_table(
            results / "coverage" / kind / f"{sample}.tsv",
            [
                {
                    heading: target,
                    "fixture Covered Fraction": "0.1",
                    "fixture Read Count": "2",
                    "fixture Mean": "0.1",
                    "fixture Trimmed Mean": "0.1",
                    "fixture Covered Bases": "10",
                    "fixture Length": "100",
                }
            ],
        )
        flagstat = results / "bam" / kind / f"{sample}.flagstat.txt"
        flagstat.parent.mkdir(parents=True)
        flagstat.write_text("2 + 0 mapped (100.00% : N/A)\n")
    (provenance / "trim").mkdir(parents=True)
    (provenance / "trim" / f"{sample}.json").write_text(
        json.dumps({"pairs_after_fastp": 100})
    )
    manifest, eligibility, samples, inventory, qc = [
        tmp_path / f"{name}.tsv"
        for name in ("manifest", "eligibility", "samples", "inventory", "qc")
    ]
    write_table(
        manifest, [{"catalog_id": "BAC00001", "role": "mag", "taxonomy": "d__Bacteria"}]
    )
    write_table(eligibility, [{"sample_id": sample, "retained_pairs": "25"}])
    write_table(samples, [{"sample_id": sample, "host_route": "none"}])
    write_table(inventory, [{"sample_id": sample, "availability": availability}])
    write_table(qc, [{"sample_id": sample, "pairs_after_fastp": "100"}])
    bacterial, viral, libraries = [
        tmp_path / f"{name}.tsv" for name in ("bacterial", "viral", "libraries")
    ]
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "workflow/scripts/aggregate_mapping.py"),
            "--samples",
            sample,
            "--results",
            str(results),
            "--provenance",
            str(provenance),
            "--manifest",
            str(manifest),
            "--eligibility",
            str(eligibility),
            "--samples-table",
            str(samples),
            "--read-inventory",
            str(inventory),
            "--bacterial",
            str(bacterial),
            "--viral",
            str(viral),
            "--libraries",
            str(libraries),
        ],
        check=True,
        capture_output=True,
    )
    with libraries.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["read_source"] == expected
    assert row["input_pairs"] == "100"  # never use the assembly subset's 25 pairs
    fasta, done = tmp_path / "viral.fna", tmp_path / "done"
    fasta.write_text(">VOTU00001\nACGT\n")
    command = [
        sys.executable,
        str(ROOT / "workflow/scripts/validate_mapping.py"),
        "--samples",
        sample,
        "--bacterial",
        str(bacterial),
        "--viral",
        str(viral),
        "--libraries",
        str(libraries),
        "--manifest",
        str(manifest),
        "--viral-fasta",
        str(fasta),
        "--library-qc",
        str(qc),
        "--output",
        str(done),
    ]
    subprocess.run(command, check=True, capture_output=True)
    assert "status\tPASS" in done.read_text()
    write_table(qc, [{"sample_id": sample, "pairs_after_fastp": "101"}])
    assert subprocess.run(command, check=False, capture_output=True).returncode != 0
