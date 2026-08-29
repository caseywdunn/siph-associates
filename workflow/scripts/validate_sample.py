import gzip
import json
import os
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, fastq_pair_count, file_record

ensure_parents([snakemake.output[0], snakemake.log[0]])
errors = []

try:
    trimmed_pairs = fastq_pair_count(snakemake.input.trim_r1, snakemake.input.trim_r2)
except Exception as exc:
    errors.append(f"trimmed FASTQ: {exc}")
    trimmed_pairs = None
try:
    host_pairs = fastq_pair_count(snakemake.input.host_r1, snakemake.input.host_r2)
except Exception as exc:
    errors.append(f"host-handled FASTQ: {exc}")
    host_pairs = None

try:
    fastp = json.loads(Path(snakemake.input.fastp).read_text())
    if int(fastp["summary"]["after_filtering"]["total_reads"]) <= 0:
        errors.append("fastp reports zero reads")
except Exception as exc:
    errors.append(f"fastp report: {exc}")

for label, path in (("Kraken2", snakemake.input.kraken), ("Bracken genus", snakemake.input.genus),
                    ("Bracken species", snakemake.input.species), ("sylph", snakemake.input.sylph)):
    if not Path(path).is_file() or Path(path).stat().st_size == 0:
        errors.append(f"{label} output is empty or absent")

try:
    with tarfile.open(snakemake.input.phyloflash, "r:gz") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
    if not members:
        errors.append("phyloFlash archive contains no files")
except Exception as exc:
    errors.append(f"phyloFlash archive: {exc}")

try:
    host_metrics = json.loads(Path(snakemake.input.host_metrics).read_text())
    if int(host_metrics["retained_pairs"]) != host_pairs:
        errors.append("host retained-pair count disagrees with metrics")
except Exception as exc:
    errors.append(f"host metrics: {exc}")

Path(snakemake.log[0]).write_text("\n".join([f"errors={len(errors)}", *errors]) + "\n")
if errors:
    raise ValueError("; ".join(errors))
atomic_json(snakemake.output[0], {
    "status": "PASS",
    "sample_id": str(snakemake.wildcards.sample),
    "study": str(snakemake.params.study),
    "host_route": str(snakemake.params.route),
    "trimmed_pairs": trimmed_pairs,
    "host_retained_pairs": host_pairs,
    "validated_outputs": [file_record(path) for path in list(snakemake.input)],
})
