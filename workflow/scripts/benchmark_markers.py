import os
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, atomic_move, ensure_parents, executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
prefix = Path(str(snakemake.params.prefix))
barrnap = executable(prefix, "barrnap")
environment = os.environ.copy()
environment["PATH"] = str(prefix / "bin") + os.pathsep + environment.get("PATH", "")
counts = Counter()

with tempfile.TemporaryDirectory(prefix=f"markers.{snakemake.wildcards.benchmark_id}.",
                                 dir=str(Path(snakemake.output.gff).parent)) as temporary:
    combined = Path(temporary) / "markers.gff"
    commands = []
    with open(combined, "w") as output, open(snakemake.log[0], "w") as log:
        for kingdom in ("bac", "arc", "euk"):
            command = [barrnap, "--kingdom", kingdom, "--threads", str(snakemake.threads), str(snakemake.input.contigs)]
            commands.append(command)
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=log, text=True, env=environment)
            if result.returncode:
                raise subprocess.CalledProcessError(result.returncode, command)
            for line in result.stdout.splitlines():
                if line and not line.startswith("#"):
                    counts[kingdom] += 1
                    fields = line.split("\t")
                    fields[-1] = fields[-1].rstrip(";") + f";kingdom={kingdom}"
                    output.write("\t".join(fields) + "\n")
    atomic_move(combined, snakemake.output.gff)

atomic_json(snakemake.output.metrics, {
    "benchmark_id": str(snakemake.wildcards.benchmark_id),
    "bacterial_rrna": counts["bac"],
    "archaeal_rrna": counts["arc"],
    "eukaryotic_rrna": counts["euk"],
    "total_rrna": sum(counts.values()),
    "software": version([barrnap, "--version"]),
    "assembly_provenance_sha256": sha256(snakemake.input.upstream),
})
