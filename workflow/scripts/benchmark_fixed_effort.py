import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, fastq_pair_count, file_record, sha256, version

if str(snakemake.params.strategy) != "fixed_effort_trimmed":
    raise ValueError(f"fixed-effort rule received strategy {snakemake.params.strategy}")
ensure_parents(list(snakemake.output) + [snakemake.log[0]])
reformat = shutil.which("reformat.sh")
if not reformat:
    raise FileNotFoundError("reformat.sh is absent from the activated environment")

with tempfile.TemporaryDirectory(prefix=f"fixed.{snakemake.wildcards.benchmark_id}.",
                                 dir=str(Path(snakemake.output.r1).parent)) as temporary:
    temporary = Path(temporary)
    r1, r2 = temporary / "R1.fastq.gz", temporary / "R2.fastq.gz"
    command = [
        reformat, f"in1={snakemake.input.r1}", f"in2={snakemake.input.r2}",
        f"out1={r1}", f"out2={r2}",
        f"samplereadstarget={int(snakemake.params.pairs)}",
        f"sampleseed={int(snakemake.params.seed)}", f"threads={snakemake.threads}",
        "overwrite=t",
    ]
    with open(snakemake.log[0], "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    pairs = fastq_pair_count(r1, r2)
    if pairs <= 0 or pairs > int(snakemake.params.pairs):
        raise ValueError(f"invalid fixed-effort pair count: {pairs}")
    atomic_move(r1, snakemake.output.r1)
    atomic_move(r2, snakemake.output.r2)

atomic_json(snakemake.output.provenance, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "strategy": "fixed_effort_trimmed",
    "target_pairs": int(snakemake.params.pairs),
    "retained_pairs": pairs,
    "seed": int(snakemake.params.seed),
    "command": command,
    "software": version([reformat, "--version"]),
    "upstream_provenance_sha256": sha256(snakemake.input.upstream),
    "outputs": [file_record(snakemake.output.r1), file_record(snakemake.output.r2)],
})
