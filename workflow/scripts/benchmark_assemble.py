import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, executable, file_record, sha256, version


def fasta_metrics(path):
    lengths, current = [], 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    total, running, n50 = sum(lengths), 0, 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= total / 2:
            n50 = length
            break
    return {"contigs": len(lengths), "total_bases": total, "max_contig": max(lengths, default=0), "n50": n50}


ensure_parents(list(snakemake.output) + [snakemake.log[0]])
prefix = Path(str(snakemake.params.prefix))
megahit = executable(prefix, "megahit")
environment = os.environ.copy()
environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")
scratch_parent = Path(snakemake.config["scratch_root"]) / "phase3_benchmark" / "megahit_tmp"
scratch_parent.mkdir(parents=True, exist_ok=True)

with tempfile.TemporaryDirectory(prefix=f"megahit.{snakemake.wildcards.benchmark_id}.",
                                 dir=str(scratch_parent)) as temporary:
    output_dir = Path(temporary) / "assembly"
    command = [megahit, "-1", str(snakemake.input.r1), "-2", str(snakemake.input.r2),
               "-t", str(snakemake.threads), "-o", str(output_dir),
               "--min-contig-len", str(snakemake.params.min_contig)]
    with open(snakemake.log[0], "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, env=environment)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    contigs = output_dir / "final.contigs.fa"
    if not contigs.is_file():
        raise FileNotFoundError("MEGAHIT did not produce final.contigs.fa")
    metrics = fasta_metrics(contigs)
    if metrics["contigs"] == 0:
        raise ValueError("MEGAHIT produced no contigs at the locked minimum length")
    atomic_move(contigs, snakemake.output.contigs)

atomic_json(snakemake.output.metrics, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "sample_id": str(snakemake.params.sample),
    "strategy": str(snakemake.params.strategy),
    "min_contig_length": int(snakemake.params.min_contig),
    **metrics,
})
atomic_json(snakemake.output.provenance, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "sample_id": str(snakemake.params.sample),
    "strategy": str(snakemake.params.strategy),
    "command": command,
    "software": version([megahit, "--version"]),
    "upstream_provenance_sha256": sha256(snakemake.input.upstream),
    "output": file_record(snakemake.output.contigs, checksum=True),
})
