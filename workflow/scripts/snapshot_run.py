import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import ensure_parents, sha256

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
pairs = [
    (snakemake.input.config, snakemake.output.config),
    (snakemake.input.samples, snakemake.output.samples),
    (snakemake.input.freeze, snakemake.output.manifest_freeze),
]
for source, destination in pairs:
    temporary = str(destination) + ".tmp"
    shutil.copyfile(source, temporary)
    os.replace(temporary, destination)
checksum_temp = str(snakemake.output.checksums) + ".tmp"
with open(checksum_temp, "w") as handle:
    for _, destination in pairs:
        handle.write(f"{sha256(destination)}  {destination}\n")
os.replace(checksum_temp, snakemake.output.checksums)
Path(snakemake.log[0]).write_text(f"snapshotted={len(pairs)}\n")
