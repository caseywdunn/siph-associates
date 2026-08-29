import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, fastq_pair_count, file_record, resolve_executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
route = str(snakemake.params.route)

with tempfile.TemporaryDirectory(prefix=f"host.{snakemake.wildcards.sample}.", dir=str(Path(snakemake.output.r1).parent)) as temp:
    temp = Path(temp)
    out1, out2 = temp / "nonhost_R1.fastq.gz", temp / "nonhost_R2.fastq.gz"
    command = []
    software = {}
    if route == "none":
        shutil.copyfile(snakemake.input.classified_r1, out1)
        shutil.copyfile(snakemake.input.classified_r2, out2)
        command = ["copy", "kraken2_classified_pairs"]
        software = {"method": "Kraken2 classified-out reuse"}
        Path(snakemake.log[0]).write_text("route=none method=kraken2_classified_out\n")
    else:
        bwa = resolve_executable(snakemake.params.prefix, "bwa")
        samtools = resolve_executable(snakemake.params.prefix, "samtools")
        command = [bwa, "mem", "-t", str(snakemake.threads), str(snakemake.params.index_prefix),
                   str(snakemake.input.r1), str(snakemake.input.r2)]
        fastq_command = [samtools, "fastq", "-f", "12", "-1", str(out1), "-2", str(out2),
                         "-0", "/dev/null", "-s", "/dev/null", "-n", "-"]
        with open(snakemake.log[0], "w") as log:
            mapper = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
            assert mapper.stdout is not None
            converter = subprocess.Popen(fastq_command, stdin=mapper.stdout, stdout=log, stderr=log)
            mapper.stdout.close()
            converter_returncode = converter.wait()
            mapper_returncode = mapper.wait()
        if mapper_returncode or converter_returncode:
            raise RuntimeError(f"host-mapping pipeline failed: bwa={mapper_returncode}, samtools={converter_returncode}")
        software = {"bwa": version([bwa]), "samtools": version([samtools, "--version"])}
        command = {"mapping": command, "both_unmapped": fastq_command}

    pairs = fastq_pair_count(out1, out2)
    os.replace(out1, snakemake.output.r1)
    os.replace(out2, snakemake.output.r2)

atomic_json(snakemake.output.metrics, {
    "sample_id": str(snakemake.wildcards.sample),
    "host_route": route,
    "retained_pairs": pairs,
    "empty_is_valid": pairs == 0,
})
atomic_json(snakemake.output.provenance, {
    "sample_id": str(snakemake.wildcards.sample),
    "host_route": route,
    "command": command,
    "software": software,
    "trim_provenance_sha256": sha256(snakemake.input.trim_provenance),
    "kraken_provenance_sha256": sha256(snakemake.input.kraken_provenance),
    "outputs": [file_record(snakemake.output.r1, checksum=True), file_record(snakemake.output.r2, checksum=True)],
})
