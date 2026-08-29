#!/usr/bin/env python3
"""Exercise Phase-2 validation and aggregation on the completed Phase-1 fixture."""
from __future__ import annotations

import csv
import json
import runpy
import tempfile
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / "tests" / "work"
SCRATCH = ROOT / "tests" / "scratch"


class NamedItems(list):
    def __init__(self, **items):
        super().__init__(items.values())
        for key, value in items.items():
            setattr(self, key, value)


def execute(script, mock):
    runpy.run_path(str(ROOT / "workflow" / "scripts" / script), init_globals={"snakemake": mock})


with (ROOT / "tests" / "fixtures" / "samples.tsv").open(newline="") as handle:
    samples = list(csv.DictReader(handle, delimiter="\t"))

with tempfile.TemporaryDirectory(prefix="phase2-fixture.") as temporary:
    temporary = Path(temporary)
    validation_paths = []
    for row in samples:
        sample = row["sample_id"]
        output = temporary / "validation" / f"{sample}.json"
        log = temporary / "logs" / f"{sample}.log"
        inputs = NamedItems(
            trim_r1=str(SCRATCH / "trimmed" / f"{sample}_R1.fastq.gz"),
            trim_r2=str(SCRATCH / "trimmed" / f"{sample}_R2.fastq.gz"),
            fastp=str(WORK / "trim" / f"{sample}.fastp.json"),
            trim_provenance=str(WORK / "provenance" / "trim" / f"{sample}.json"),
            kraken=str(WORK / "screens" / "kraken" / f"{sample}.report.tsv"),
            genus=str(WORK / "screens" / "bracken" / f"{sample}.G.tsv"),
            species=str(WORK / "screens" / "bracken" / f"{sample}.S.tsv"),
            sylph=str(WORK / "screens" / "sylph" / f"{sample}.profile.tsv"),
            phyloflash=str(WORK / "screens" / "phyloflash" / f"{sample}.tar.gz"),
            kraken_provenance=str(WORK / "provenance" / "kraken_bracken" / f"{sample}.json"),
            sylph_provenance=str(WORK / "provenance" / "sylph" / f"{sample}.json"),
            phyloflash_provenance=str(WORK / "provenance" / "phyloflash" / f"{sample}.json"),
        )
        mock = SimpleNamespace(
            config={"manifest": str(ROOT / "manifest.csv")}, input=inputs, output=[str(output)],
            log=[str(log)], wildcards=SimpleNamespace(sample=sample),
            params=SimpleNamespace(
                library_id=row["library_id"], specimen_id=row["specimen_id"], study=row["study"],
                route=row["host_route"], raw_pairs=50000, cap_pairs=400000000,
            ),
        )
        execute("validate_screen_sample.py", mock)
        validation_paths.append(str(output))

    ids = [row["sample_id"] for row in samples]
    aggregate_inputs = SimpleNamespace(
        validations=validation_paths,
        genus=[str(WORK / "screens" / "bracken" / f"{sample}.G.tsv") for sample in ids],
        species=[str(WORK / "screens" / "bracken" / f"{sample}.S.tsv") for sample in ids],
        sylph=[str(WORK / "screens" / "sylph" / f"{sample}.profile.tsv") for sample in ids],
        phyloflash=[str(WORK / "screens" / "phyloflash" / f"{sample}.tar.gz") for sample in ids],
    )
    qc = temporary / "aggregation" / "library_qc.tsv"
    nominations = temporary / "aggregation" / "candidate_nominations.tsv"
    provenance = temporary / "aggregation" / "provenance.json"
    aggregate_output = NamedItems(qc=str(qc), nominations=str(nominations), provenance=str(provenance))
    aggregate_mock = SimpleNamespace(
        config={"manifest": str(ROOT / "manifest.csv")}, input=aggregate_inputs,
        output=aggregate_output, log=[str(temporary / "logs" / "aggregate.log")],
        params=SimpleNamespace(sample_ids=ids), wildcards=SimpleNamespace(scope="fixture"),
    )
    execute("aggregate_phase2.py", aggregate_mock)

    stage = temporary / "screen_fixture.done"
    stage_inputs = SimpleNamespace(
        validations=validation_paths, qc=str(qc), nominations=str(nominations), provenance=str(provenance),
    )
    stage_mock = SimpleNamespace(
        config={"manifest": str(ROOT / "manifest.csv")}, input=stage_inputs,
        output=[str(stage)], log=[str(temporary / "logs" / "stage.log")],
        params=SimpleNamespace(expected_samples=3, expected_studies=3, expected_routes=3),
    )
    execute("validate_phase2_stage.py", stage_mock)
    result = dict(line.split("\t", 1) for line in stage.read_text().splitlines())
    if result.get("status") != "PASS" or result.get("samples") != "3":
        raise AssertionError(result)
    if sum(1 for _ in qc.open()) != 4:
        raise AssertionError("fixture QC aggregate does not have header plus three rows")
    print(stage.read_text(), end="")
