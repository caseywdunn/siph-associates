import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, file_record, resolve_executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
kraken = resolve_executable(snakemake.params.prefix, "kraken2")
bracken = resolve_executable(snakemake.params.prefix, "bracken")
fastp_report = json.loads(Path(snakemake.input.fastp).read_text())
read_length = round(float(fastp_report["summary"]["after_filtering"]["read1_mean_length"]))
database = str(snakemake.input.database)

with tempfile.TemporaryDirectory(prefix=f"kraken.{snakemake.wildcards.sample}.", dir=str(Path(snakemake.output.report).parent)) as temp:
    temp = Path(temp)
    report = temp / "report.tsv"
    assignments = temp / "assignments.tsv.gz"
    classified_template = str(temp / "classified_#.fastq")
    command = [
        kraken, "--db", database, "--threads", str(snakemake.threads), "--paired",
        "--gzip-compressed", "--report", str(report), "--report-minimizer-data",
        "--classified-out", classified_template, str(snakemake.input.r1), str(snakemake.input.r2),
    ]
    with open(snakemake.log[0], "w") as log, gzip.open(assignments, "wb") as output:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
        assert process.stdout is not None
        shutil.copyfileobj(process.stdout, output)
        process.stdout.close()
        returncode = process.wait()
    if returncode:
        raise subprocess.CalledProcessError(returncode, command)
    if not report.is_file() or report.stat().st_size == 0:
        raise ValueError("Kraken2 report is empty")

    classified = []
    for mate in (1, 2):
        plain = temp / f"classified_{mate}.fastq"
        compressed = temp / f"classified_{mate}.fastq.gz"
        with open(plain, "rb") as source, gzip.open(compressed, "wb") as target:
            shutil.copyfileobj(source, target)
        classified.append(compressed)

    bracken_outputs = {}
    with open(snakemake.log[0], "a") as log:
        for level, final in (("G", snakemake.output.genus), ("S", snakemake.output.species)):
            output = temp / f"bracken.{level}.tsv"
            kreport = temp / f"bracken.{level}.kreport"
            cmd = [bracken, "-d", database, "-i", str(report), "-o", str(output),
                   "-w", str(kreport), "-r", str(read_length), "-l", level,
                   "-t", str(snakemake.params.threshold)]
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
            if result.returncode:
                raise subprocess.CalledProcessError(result.returncode, cmd)
            if not output.is_file() or output.stat().st_size == 0:
                raise ValueError(f"Bracken {level} output is empty")
            bracken_outputs[final] = output

    moves = {
        snakemake.output.report: report,
        snakemake.output.assignments: assignments,
        snakemake.output.classified_r1: classified[0],
        snakemake.output.classified_r2: classified[1],
        **bracken_outputs,
    }
    for final, temporary in moves.items():
        os.replace(temporary, final)

atomic_json(snakemake.output.provenance, {
    "sample_id": str(snakemake.wildcards.sample),
    "command": command,
    "software": {"kraken2": version([kraken, "--version"]), "bracken": version([bracken, "-v"])},
    "parameters": {"read_length": read_length, "bracken_threshold": int(snakemake.params.threshold)},
    "database": file_record(Path(database).resolve() / "hash.k2d"),
    "upstream_provenance_sha256": sha256(snakemake.input.trim_provenance),
    "outputs": [file_record(path) for path in list(snakemake.output)[:-1]],
})
