import os
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, file_record, resolve_executable, sha256, version

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
phyloflash = resolve_executable(snakemake.params.prefix, "phyloFlash.pl")
library = str(snakemake.wildcards.sample)

with tempfile.TemporaryDirectory(prefix=f"phyloflash.{library}.", dir=str(Path(snakemake.output.archive).parent)) as temp:
    temp = Path(temp)
    work = temp / "work"
    work.mkdir()
    archive = temp / "output.tar.gz"
    command = [phyloflash, "-lib", library, "-read1", str(snakemake.input.r1),
               "-read2", str(snakemake.input.r2), "-dbhome", str(snakemake.input.database),
               "-CPUs", str(snakemake.threads), "-almosteverything", "-log"]
    environment = os.environ.copy()
    prefix = str(snakemake.params.prefix)
    environment["PATH"] = str(Path(prefix) / "bin") + os.pathsep + environment.get("PATH", "")
    with open(snakemake.log[0], "w") as log:
        result = subprocess.run(command, cwd=work, env=environment, stdout=log,
                                stderr=subprocess.STDOUT, text=True)
    if result.returncode:
        raise subprocess.CalledProcessError(result.returncode, command)
    files = [path for path in work.rglob("*") if path.is_file()]
    if not files:
        raise ValueError("phyloFlash completed without producing files")
    with tarfile.open(archive, "w:gz") as tar:
        for path in sorted(files):
            tar.add(path, arcname=path.relative_to(work))
    os.replace(archive, snakemake.output.archive)

atomic_json(snakemake.output.provenance, {
    "sample_id": library,
    "command": command,
    "software": version([phyloflash, "-version"]),
    "database": str(snakemake.input.database),
    "upstream_provenance_sha256": sha256(snakemake.input.trim_provenance),
    "archive": file_record(snakemake.output.archive, checksum=True),
})
