import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, fastq_pair_count, file_record, executable, sha256, version

if str(snakemake.params.strategy) != "reference_depleted":
    raise ValueError(f"reference-depletion rule received strategy {snakemake.params.strategy}")
if str(snakemake.params.route) == "none":
    raise ValueError("reference depletion cannot run for host_route=none")
ensure_parents(list(snakemake.output) + [snakemake.log[0]])
prefix = Path(str(snakemake.params.prefix))
bwa, samtools = executable(prefix, "bwa"), executable(prefix, "samtools")
environment = os.environ.copy()
environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")

with tempfile.TemporaryDirectory(prefix=f"refdep.{snakemake.wildcards.benchmark_id}.",
                                 dir=str(Path(snakemake.output.r1).parent)) as temporary:
    temporary = Path(temporary)
    r1, r2 = temporary / "R1.fastq.gz", temporary / "R2.fastq.gz"
    mapping = [bwa, "mem", "-t", str(snakemake.threads), str(snakemake.params.index_prefix),
               str(snakemake.input.r1), str(snakemake.input.r2)]
    filtering = [samtools, "fastq", "-f", "12", "-1", str(r1), "-2", str(r2),
                 "-0", "/dev/null", "-s", "/dev/null", "-n", "-"]
    with open(snakemake.log[0], "w") as log:
        mapper = subprocess.Popen(mapping, stdout=subprocess.PIPE, stderr=log, env=environment)
        assert mapper.stdout is not None
        converter = subprocess.Popen(filtering, stdin=mapper.stdout, stdout=log, stderr=log, env=environment)
        mapper.stdout.close()
        converter_code = converter.wait()
        mapper_code = mapper.wait()
    if mapper_code or converter_code:
        raise RuntimeError(f"reference-depletion pipeline failed: bwa={mapper_code}, samtools={converter_code}")
    pairs = fastq_pair_count(r1, r2)
    atomic_move(r1, snakemake.output.r1)
    atomic_move(r2, snakemake.output.r2)

atomic_json(snakemake.output.metrics, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "sample_id": str(snakemake.params.sample),
    "host_route": str(snakemake.params.route),
    "retained_pairs": pairs,
    "empty_is_valid": pairs == 0,
})
atomic_json(snakemake.output.provenance, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "strategy": "reference_depleted",
    "command": {"mapping": mapping, "both_unmapped": filtering},
    "software": {"bwa": version([bwa]), "samtools": version([samtools, "--version"])},
    "upstream_provenance_sha256": sha256(snakemake.input.upstream),
    "outputs": [file_record(snakemake.output.r1), file_record(snakemake.output.r2)],
})
