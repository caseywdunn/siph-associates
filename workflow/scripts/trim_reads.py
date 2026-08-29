import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, file_record, sha256, version


def concatenate(inputs, output):
    with open(output, "wb") as target:
        for path in inputs:
            with open(path, "rb") as source:
                shutil.copyfileobj(source, target, 1024 * 1024)


def deterministic_cap(r1_paths, r2_paths, keep, seed, out1, out2, threads, log):
    """Create an exact seeded paired subsample with BBTools."""
    reformat = shutil.which("reformat.sh")
    if not reformat:
        raise FileNotFoundError("deterministic capping requires reformat.sh from the declared BBMap module")
    raw1, raw2 = out1.parent / "uncapped_R1.fastq.gz", out1.parent / "uncapped_R2.fastq.gz"
    concatenate(r1_paths, raw1)
    concatenate(r2_paths, raw2)
    command = [
        reformat, f"in1={raw1}", f"in2={raw2}", f"out1={out1}", f"out2={out2}",
        f"samplereadstarget={keep}", f"sampleseed={seed}",
        f"threads={max(1, threads // 2)}", "overwrite=t",
    ]
    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    raw1.unlink()
    raw2.unlink()
    return command


ensure_parents(list(snakemake.output) + [snakemake.log[0]])
expected = int(snakemake.params.expected_pairs)
cap = int(snakemake.params.cap_pairs)
kept = min(expected, cap)
fastp = shutil.which("fastp")
if not fastp:
    raise FileNotFoundError("fastp is absent from PATH; run Snakemake with --use-envmodules")

with tempfile.TemporaryDirectory(prefix=f"trim.{snakemake.wildcards.sample}.", dir=str(Path(snakemake.output.r1).parent)) as temp:
    temp = Path(temp)
    staged1, staged2 = temp / "input_R1.fastq.gz", temp / "input_R2.fastq.gz"
    out1, out2 = temp / "trimmed_R1.fastq.gz", temp / "trimmed_R2.fastq.gz"
    report_json, report_html = temp / "fastp.json", temp / "fastp.html"
    with open(snakemake.log[0], "w") as log:
        if expected > cap:
            cap_command = deterministic_cap(
                snakemake.input.r1, snakemake.input.r2, cap, int(snakemake.params.cap_seed),
                staged1, staged2, int(snakemake.threads), log,
            )
            cap_method = "bbtools_exact_seeded_pair_sampling"
        else:
            concatenate(snakemake.input.r1, staged1)
            concatenate(snakemake.input.r2, staged2)
            cap_command = None
            cap_method = "not_capped"
        command = [
            fastp, "-i", str(staged1), "-I", str(staged2), "-o", str(out1), "-O", str(out2),
            "--thread", str(snakemake.threads), "--detect_adapter_for_pe", "--dont_eval_duplication",
            "--json", str(report_json), "--html", str(report_html),
        ]
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    report = json.loads(report_json.read_text())
    passed_reads = int(report["summary"]["after_filtering"]["total_reads"])
    if passed_reads <= 0 or passed_reads % 2:
        raise ValueError(f"fastp produced invalid paired read count: {passed_reads}")
    for temporary, final in ((out1, snakemake.output.r1), (out2, snakemake.output.r2),
                             (report_json, snakemake.output.fastp_json), (report_html, snakemake.output.fastp_html)):
        atomic_move(temporary, final)

atomic_json(snakemake.output.provenance, {
    "sample_id": str(snakemake.wildcards.sample),
    "command": {"cap": cap_command, "fastp": command},
    "software": version([fastp, "--version"]),
    "parameters": {"cap_pairs": cap, "cap_seed": int(snakemake.params.cap_seed), "cap_method": cap_method},
    "run_snapshot_sha256": sha256(snakemake.input.run_snapshot),
    "raw_pairs_expected": expected,
    "pairs_entering_fastp": kept,
    "pairs_after_fastp": passed_reads // 2,
    "inputs": [file_record(path, checksum=True) for path in list(snakemake.input.r1) + list(snakemake.input.r2)],
    "outputs": [file_record(snakemake.output.r1, checksum=True), file_record(snakemake.output.r2, checksum=True)],
})
