import csv
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, executable, sha256, version


def fasta_metrics(path):
    records, bases = 0, 0
    if path.is_file():
        with path.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    records += 1
                else:
                    bases += len(line.strip())
    return records, bases


ensure_parents(list(snakemake.output) + [snakemake.log[0]])
prefix = Path(str(snakemake.params.prefix))
genomad = executable(prefix, "genomad")
environment = os.environ.copy()
environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")

with tempfile.TemporaryDirectory(prefix=f"genomad.{snakemake.wildcards.benchmark_id}.",
                                 dir=str(Path(snakemake.output.viruses).parent)) as temporary:
    temporary = Path(temporary)
    input_fasta = temporary / "assembly.fasta"
    os.link(snakemake.input.contigs, input_fasta)
    output_dir = temporary / "genomad"
    command = [genomad, "end-to-end", "--threads", str(snakemake.threads), "--cleanup",
               str(input_fasta), str(output_dir), str(snakemake.input.database)]
    with open(snakemake.log[0], "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True, env=environment)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    summary_dir = output_dir / "assembly_summary"
    virus_fasta = summary_dir / "assembly_virus.fna"
    virus_summary = summary_dir / "assembly_virus_summary.tsv"
    if not virus_fasta.is_file():
        virus_fasta.write_text("")
    if not virus_summary.is_file():
        virus_summary.write_text("seq_name\n")
    records, bases = fasta_metrics(virus_fasta)
    atomic_move(virus_fasta, snakemake.output.viruses)
    atomic_move(virus_summary, snakemake.output.summary)

atomic_json(snakemake.output.metrics, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "viral_contigs": records,
    "viral_bases": bases,
    "software": version([genomad, "--version"]),
    "assembly_provenance_sha256": sha256(snakemake.input.upstream),
})
