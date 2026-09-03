import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, executable, fastq_pair_count, file_record, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
strategy = str(snakemake.params.strategy)
route = str(snakemake.params.route)
prefix = Path(str(snakemake.params.prefix))

with tempfile.TemporaryDirectory(prefix=f"input.{snakemake.wildcards.sample}.",
                                 dir=str(Path(snakemake.output.r1).parent)) as temporary:
    temporary = Path(temporary)
    r1, r2 = temporary / "R1.fastq.gz", temporary / "R2.fastq.gz"
    if strategy == "fixed_effort_trimmed":
        if route != "none":
            raise ValueError(f"fixed-effort strategy assigned to reference route {route}")
        reformat = shutil.which("reformat.sh")
        if not reformat:
            raise FileNotFoundError("reformat.sh is absent from the activated environment")
        command = [
            reformat, f"in1={snakemake.input.r1}", f"in2={snakemake.input.r2}",
            f"out1={r1}", f"out2={r2}",
            f"samplereadstarget={int(snakemake.params.fixed_pairs)}",
            f"sampleseed={int(snakemake.params.fixed_seed)}", f"threads={snakemake.threads}",
            "overwrite=t",
        ]
        with open(snakemake.log[0], "w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
        if result.returncode:
            raise subprocess.CalledProcessError(result.returncode, command)
        software = {"bbmap": version([reformat, "--version"])}
    elif strategy == "reference_depleted":
        if route == "none":
            raise ValueError("reference-depleted strategy assigned without a host reference")
        bwa, samtools = executable(prefix, "bwa"), executable(prefix, "samtools")
        mapping = [bwa, "mem", "-t", str(snakemake.threads), str(snakemake.params.index_prefix),
                   str(snakemake.input.r1), str(snakemake.input.r2)]
        fastq = [samtools, "fastq", "-f", "12", "-1", str(r1), "-2", str(r2),
                 "-0", "/dev/null", "-s", "/dev/null", "-n", "-"]
        with open(snakemake.log[0], "w") as log:
            mapper = subprocess.Popen(mapping, stdout=subprocess.PIPE, stderr=log)
            assert mapper.stdout is not None
            converter = subprocess.Popen(fastq, stdin=mapper.stdout, stdout=log, stderr=log)
            mapper.stdout.close()
            converter_code = converter.wait()
            mapper_code = mapper.wait()
        if mapper_code or converter_code:
            raise RuntimeError(f"host depletion failed: bwa={mapper_code}, samtools={converter_code}")
        command = {"mapping": mapping, "both_unmapped": fastq}
        software = {"bwa": version([bwa]), "samtools": version([samtools, "--version"])}
    else:
        raise ValueError(f"unknown Phase-3 cohort strategy: {strategy}")

    retained_pairs = fastq_pair_count(r1, r2)
    if retained_pairs <= 0:
        raise ValueError("Phase-3 cohort input contains no read pairs")
    os.replace(r1, snakemake.output.r1)
    os.replace(r2, snakemake.output.r2)

upstream = json.loads(Path(snakemake.input.upstream).read_text())
input_pairs = int(upstream["pairs_after_fastp"])
atomic_json(snakemake.output.metrics, {
    "sample_id": str(snakemake.wildcards.sample),
    "host_route": route,
    "strategy": strategy,
    "input_pairs": input_pairs,
    "retained_pairs": retained_pairs,
    "retained_fraction": retained_pairs / input_pairs,
})
atomic_json(snakemake.output.provenance, {
    "sample_id": str(snakemake.wildcards.sample),
    "host_route": route,
    "strategy": strategy,
    "fixed_effort_pairs": int(snakemake.params.fixed_pairs) if strategy == "fixed_effort_trimmed" else None,
    "fixed_effort_seed": int(snakemake.params.fixed_seed) if strategy == "fixed_effort_trimmed" else None,
    "command": command,
    "software": software,
    "upstream_provenance_sha256": sha256(snakemake.input.upstream),
    "outputs": [file_record(snakemake.output.r1), file_record(snakemake.output.r2)],
})
