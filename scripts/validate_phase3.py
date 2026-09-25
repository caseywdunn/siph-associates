#!/usr/bin/env python3
"""Validate the frozen Phase-3 benchmark design and optional completed outputs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REPORT = ROOT / "data" / "metadata" / "phase3_validation.txt"
parser = argparse.ArgumentParser()
parser.add_argument("--scope", choices=("static", "benchmark", "assembly"), default="static")
args = parser.parse_args()
errors = []


def require(condition, message):
    if not condition:
        errors.append(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


samples = {row["sample_id"]: row for row in table(ROOT / "config" / "samples.tsv")
           if row["include_primary"].lower() == "true"}
panel = table(ROOT / "config" / "phase3_benchmark.tsv")
settings = json.loads((ROOT / "config" / "phase3_benchmark.json").read_text())
cohort_settings = json.loads((ROOT / "config" / "phase3_cohort.json").read_text())
assembly_settings = json.loads((ROOT / "config" / "phase3_assembly.json").read_text())
ids = [row["benchmark_id"] for row in panel]
panel_samples = sorted({row["sample_id"] for row in panel})
by_sample = defaultdict(set)
for row in panel:
    by_sample[row["sample_id"]].add(row["strategy"])

require(len(panel) == 12, f"expected 12 strategy rows, found {len(panel)}")
require(len(ids) == len(set(ids)), "benchmark IDs are duplicated")
require(len(panel_samples) == 6, f"expected six benchmark samples, found {len(panel_samples)}")
require(set(panel_samples).issubset(samples), "benchmark contains an unknown sample")
require({row["signal_stratum"] for row in panel} == {"high", "low"}, "benchmark lacks high/low signal strata")
require({samples[s]["study"] for s in panel_samples} == {"Church2025", "Ahuja2024", "Ahuja2026"},
        "benchmark does not cover all studies")
require({samples[s]["host_route"] for s in panel_samples} == {"P_physalis", "N_septata", "none"},
        "benchmark does not cover all host routes")
for sample, strategies in by_sample.items():
    route = samples[sample]["host_route"]
    expected = ({"reference_depleted", "kraken_nominated"} if route != "none"
                else {"kraken_nominated", "fixed_effort_trimmed"})
    require(strategies == expected, f"strategy pair mismatch for {sample}: {sorted(strategies)}")
require(int(settings.get("fixed_effort_pairs", 0)) == 25_000_000, "fixed-effort read ceiling is not 25M pairs")
require(int(settings.get("fixed_effort_seed", 0)) > 0, "fixed-effort seed is absent")
require(int(settings.get("min_contig_length", 0)) == 1000, "MEGAHIT minimum contig length is not 1 kb")
require(cohort_settings.get("strategy_by_host_route") == {
    "P_physalis": "reference_depleted", "N_septata": "reference_depleted",
    "none": "fixed_effort_trimmed",
}, "Phase-3 cohort route strategies do not match the accepted benchmark decision")
require(int(cohort_settings.get("fixed_effort_pairs", 0)) == 25_000_000,
        "Phase-3 cohort fixed-effort input is not 25M pairs")
require(int(cohort_settings.get("minimum_assembly_input_pairs", 0)) == 175_000,
        "Phase-3 cohort assembly eligibility floor is not 175,000 pairs")
benchmark_sentinel = ROOT / cohort_settings.get("benchmark_sentinel", "")
require(benchmark_sentinel.is_file(), "accepted Phase-3 benchmark sentinel is absent")
if benchmark_sentinel.is_file():
    benchmark_values = dict(
        line.split("\t", 1) for line in benchmark_sentinel.read_text().splitlines() if "\t" in line
    )
    require(benchmark_values.get("status") == "PASS", "Phase-3 benchmark sentinel is not PASS")

for label, path in settings.get("databases", {}).items():
    require(Path(path).exists(), f"missing {label} database: {path}")
for label, prefix in settings.get("tool_prefixes", {}).items():
    require(Path(prefix).is_dir(), f"missing {label} tool prefix: {prefix}")
required_tools = {"assembly": ("megahit", "bwa", "samtools", "metabat2", "jgi_summarize_bam_contig_depths"),
                  "barrnap": ("barrnap",), "genomad": ("genomad",), "checkm2": ("checkm2",)}
for group, tools in required_tools.items():
    prefix = Path(settings["tool_prefixes"][group]) / "bin"
    for tool in tools:
        require((prefix / tool).is_file(), f"missing executable: {prefix / tool}")

phase2 = ROOT / "data" / "results" / "stages" / "screen_cohort.done"
require(phase2.is_file(), "Phase-2 cohort sentinel is absent")
if phase2.is_file():
    values = dict(line.split("\t", 1) for line in phase2.read_text().splitlines() if "\t" in line)
    require(values.get("status") == "PASS" and values.get("samples") == "205", "Phase-2 cohort sentinel is not accepted")

for key in ("min_contig_length", "metabat_min_contig", "metabat_seed"):
    require(assembly_settings.get(key) == settings.get(key),
            f"cohort assembly {key} differs from the accepted benchmark")
for label, path in assembly_settings.get("databases", {}).items():
    require(Path(path).exists(), f"missing {label} database: {path}")
    if label in settings.get("databases", {}):
        require(path == settings["databases"][label], f"cohort {label} database differs from the benchmark")
require(assembly_settings.get("tool_prefixes") == settings.get("tool_prefixes"),
        "cohort assembly tool prefixes differ from the accepted benchmark")
require((Path(assembly_settings["tool_prefixes"]["genomad"]) / "bin" / "checkv").is_file(),
        "missing CheckV executable in the geNomad environment")
for label, resource in assembly_settings.get("resources", {}).items():
    require(resource.get("partition") != "week" or int(resource.get("runtime", 0)) >= 2880,
            f"{label} requests week for less than 48 hours")
membership_rows = table(ROOT / assembly_settings["membership"])
eligibility_path = ROOT / "data" / "results" / "phase3_cohort" / "input_eligibility.tsv"
require(eligibility_path.is_file(), "accepted Phase-3 input-eligibility table is absent")
if eligibility_path.is_file():
    eligibility = {row["sample_id"]: row for row in table(eligibility_path)}
    require({row["sample_id"] for row in membership_rows} == set(eligibility) == set(samples),
            "frozen assembly membership does not match the eligibility table and cohort")
    for row in membership_rows:
        source = eligibility.get(row["sample_id"], {})
        require(all(row[key] == source.get(key) for key in (
            "host_route", "strategy", "retained_pairs", "minimum_assembly_input_pairs", "assembly_eligible")),
            f"frozen assembly membership differs from the eligibility table: {row['sample_id']}")
assembly_ids = sorted(row["sample_id"] for row in membership_rows if row["assembly_eligible"] == "true")

snake_text = "\n".join(path.read_text() for path in [ROOT / "Snakefile", *sorted((ROOT / "workflow" / "rules").glob("*.smk"))])
for rule in ("benchmark_trim_reads", "benchmark_assembly", "benchmark_markers", "benchmark_viruses",
             "benchmark_binning", "benchmark_checkm2", "benchmark_host_carryover",
             "aggregate_phase3_benchmark", "phase3_benchmark", "phase3_prepare_cohort_input",
             "aggregate_phase3_inputs", "phase3_inputs", "assemble_contigs_megahit", "find_rrna_barrnap",
             "identify_viruses_genomad", "assess_viruses_checkv", "bin_contigs_metabat2",
             "assess_bins_checkm2", "summarize_assemblies", "validate_assemblies", "phase3_assembly"):
    require(f"rule {rule}:" in snake_text, f"missing Phase-3 rule: {rule}")

validated_rows = 0
if args.scope == "benchmark":
    decision = ROOT / "data" / "results" / "phase3_benchmark" / "decision_table.tsv"
    sentinel = ROOT / "data" / "results" / "stages" / "phase3_benchmark.done"
    require(decision.is_file(), "benchmark decision table is absent")
    require(sentinel.is_file(), "benchmark sentinel is absent")
    if decision.is_file():
        rows = table(decision)
        validated_rows = len(rows)
        require(len(rows) == len(panel), "benchmark decision-table row count mismatch")
        require({row.get("benchmark_id") for row in rows} == set(ids), "benchmark decision-table ID mismatch")
    if sentinel.is_file():
        values = dict(line.split("\t", 1) for line in sentinel.read_text().splitlines() if "\t" in line)
        require(values.get("status") == "PASS", "Phase-3 benchmark sentinel is not PASS")

validated_assemblies = 0
if args.scope == "assembly":
    results = ROOT / "data" / "results"
    summary_path = results / "phase3_cohort" / "assembly_summary.tsv"
    sentinel = results / "stages" / "phase3_assembly.done"
    require(summary_path.is_file(), "assembly summary is absent")
    require(sentinel.is_file(), "assembly sentinel is absent")
    if sentinel.is_file():
        values = dict(line.split("\t", 1) for line in sentinel.read_text().splitlines() if "\t" in line)
        require(values.get("status") == "PASS", "Phase-3 assembly sentinel is not PASS")
    if summary_path.is_file():
        summary = {row["sample_id"]: row for row in table(summary_path)}
        require(set(summary) == set(samples), "assembly summary does not cover exactly the cohort")
        for sample in samples:
            row = summary.get(sample, {})
            if sample not in assembly_ids:
                require(row.get("assembly_status") == "excluded_below_floor",
                        f"below-floor library is not a documented exclusion: {sample}")
                continue
            require(row.get("assembly_status") in ("assembled", "no_contigs"), f"library not assembled: {sample}")
            outputs = [results / "phase3_cohort" / kind / f"{sample}{suffix}" for kind, suffix in (
                ("assemblies", ".fasta"), ("markers", ".gff"), ("viruses", ".fna"), ("checkv", ".tsv"),
                ("bins", ".tar.gz"), ("bins", ".depth.tsv"), ("checkm2", ".tsv"), ("checkm2", ".status"))]
            for path in outputs:
                require(path.is_file(), f"missing assembly-branch output: {path.relative_to(ROOT)}")
            if outputs[0].is_file():
                require(sha256(outputs[0]) == row.get("assembly_sha256"),
                        f"assembly differs from its summarized checksum: {sample}")
                with outputs[0].open() as handle:
                    header = handle.readline()
                require(not header or header.startswith(f">{sample}__"),
                        f"contig IDs do not carry the source library: {sample}")
            validated_assemblies += 1

metrics = {
    "status": "PASS" if not errors else "FAIL", "scope": args.scope,
    "benchmark_rows": len(panel), "benchmark_samples": len(panel_samples),
    "validated_rows": validated_rows, "assembly_samples": len(assembly_ids),
    "validated_assemblies": validated_assemblies, "errors": len(errors),
}
REPORT.parent.mkdir(parents=True, exist_ok=True)
REPORT.write_text("".join(f"{key}\t{value}\n" for key, value in metrics.items()) +
                  "".join(f"ERROR\t{error}\n" for error in errors))
print(REPORT.read_text(), end="")
if errors:
    raise SystemExit(1)
