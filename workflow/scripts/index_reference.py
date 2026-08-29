import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, file_record, resolve_executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
bwa = resolve_executable(snakemake.params.prefix, "bwa")
final_prefix = Path(snakemake.params.index_prefix)
with tempfile.TemporaryDirectory(prefix="bwa_index.", dir=str(final_prefix.parent)) as temp:
    temp_prefix = Path(temp) / "host"
    command = [bwa, "index", "-p", str(temp_prefix), str(snakemake.input.reference)]
    with open(snakemake.log[0], "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        source = Path(str(temp_prefix) + suffix)
        if not source.is_file() or source.stat().st_size == 0:
            raise ValueError(f"BWA failed to produce nonempty {source}")
        os.replace(source, str(final_prefix) + suffix)

atomic_json(snakemake.output.provenance, {
    "route": str(snakemake.wildcards.route),
    "software": version([bwa]),
    "command": command,
    "reference": file_record(snakemake.input.reference, checksum=True),
    "run_snapshot_sha256": sha256(snakemake.input.run_snapshot),
    "outputs": [file_record(path) for path in list(snakemake.output)[:-1]],
})
