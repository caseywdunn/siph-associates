import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, file_record, sha256, version


def concatenate(inputs, output):
    with open(output, "wb") as target:
        for path in inputs:
            with open(path, "rb") as source:
                shutil.copyfileobj(source, target, 1024 * 1024)


def record(handle):
    lines = [handle.readline() for _ in range(4)]
    if not lines[0]:
        return None
    if not all(lines):
        raise ValueError("incomplete FASTQ record during deterministic capping")
    return lines


def systematic_cap(r1_paths, r2_paths, total, keep, seed, out1, out2):
    selected = 0
    seen = 0
    offset = random.Random(seed).random()
    next_index = int(offset * total / keep)
    with gzip.open(out1, "wt") as target1, gzip.open(out2, "wt") as target2:
        for path1, path2 in zip(r1_paths, r2_paths):
            with gzip.open(path1, "rt") as source1, gzip.open(path2, "rt") as source2:
                while True:
                    pair1, pair2 = record(source1), record(source2)
                    if pair1 is None and pair2 is None:
                        break
                    if pair1 is None or pair2 is None:
                        raise ValueError("mate files contain different record counts")
                    if seen == next_index and selected < keep:
                        target1.writelines(pair1)
                        target2.writelines(pair2)
                        selected += 1
                        next_index = int((selected + offset) * total / keep)
                    seen += 1
    if seen != total or selected != keep:
        raise ValueError(f"cap reconciliation failed: observed={seen}, expected={total}, selected={selected}, keep={keep}")


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
    if expected > cap:
        systematic_cap(snakemake.input.r1, snakemake.input.r2, expected, cap,
                       int(snakemake.params.cap_seed), staged1, staged2)
        cap_method = "seeded_deterministic_systematic_even_spacing"
    else:
        concatenate(snakemake.input.r1, staged1)
        concatenate(snakemake.input.r2, staged2)
        cap_method = "not_capped"
    out1, out2 = temp / "trimmed_R1.fastq.gz", temp / "trimmed_R2.fastq.gz"
    report_json, report_html = temp / "fastp.json", temp / "fastp.html"
    command = [
        fastp, "-i", str(staged1), "-I", str(staged2), "-o", str(out1), "-O", str(out2),
        "--thread", str(snakemake.threads), "--detect_adapter_for_pe", "--dont_eval_duplication",
        "--json", str(report_json), "--html", str(report_html),
    ]
    with open(snakemake.log[0], "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    report = json.loads(report_json.read_text())
    passed_reads = int(report["summary"]["after_filtering"]["total_reads"])
    if passed_reads <= 0 or passed_reads % 2:
        raise ValueError(f"fastp produced invalid paired read count: {passed_reads}")
    for temporary, final in ((out1, snakemake.output.r1), (out2, snakemake.output.r2),
                             (report_json, snakemake.output.fastp_json), (report_html, snakemake.output.fastp_html)):
        os.replace(temporary, final)

atomic_json(snakemake.output.provenance, {
    "sample_id": str(snakemake.wildcards.sample),
    "command": command,
    "software": version([fastp, "--version"]),
    "parameters": {"cap_pairs": cap, "cap_seed": int(snakemake.params.cap_seed), "cap_method": cap_method},
    "run_snapshot_sha256": sha256(snakemake.input.run_snapshot),
    "raw_pairs_expected": expected,
    "pairs_entering_fastp": kept,
    "pairs_after_fastp": passed_reads // 2,
    "inputs": [file_record(path, checksum=True) for path in list(snakemake.input.r1) + list(snakemake.input.r2)],
    "outputs": [file_record(snakemake.output.r1, checksum=True), file_record(snakemake.output.r2, checksum=True)],
})
