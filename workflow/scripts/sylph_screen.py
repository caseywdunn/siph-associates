import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, file_record, resolve_executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
sylph = resolve_executable(snakemake.params.prefix, "sylph")
databases = list(snakemake.input.databases)

with tempfile.TemporaryDirectory(prefix=f"sylph.{snakemake.wildcards.sample}.", dir=str(Path(snakemake.output.profile).parent)) as temp:
    profile = Path(temp) / "profile.tsv"
    command = [sylph, "profile", "-1", str(snakemake.input.r1), "-2", str(snakemake.input.r2),
               "-t", str(snakemake.threads), "--estimate-unknown", "-o", str(profile), *databases]
    with open(snakemake.log[0], "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    if not profile.is_file() or profile.stat().st_size == 0:
        raise ValueError("sylph profile is empty")
    os.replace(profile, snakemake.output.profile)

atomic_json(snakemake.output.provenance, {
    "sample_id": str(snakemake.wildcards.sample),
    "command": command,
    "software": version([sylph, "--version"]),
    "databases": [file_record(path) for path in databases],
    "upstream_provenance_sha256": sha256(snakemake.input.trim_provenance),
    "output": file_record(snakemake.output.profile),
})
