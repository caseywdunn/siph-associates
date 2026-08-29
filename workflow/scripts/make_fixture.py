import gzip
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, normalize_read_id, sha256

if not snakemake.params.enabled:
    raise ValueError("fixture rule requested while fixture_enabled is false")

outputs = [snakemake.output.r1, snakemake.output.r2, snakemake.output.metadata, snakemake.log[0]]
ensure_parents(outputs)
tmp1 = str(snakemake.output.r1) + ".tmp"
tmp2 = str(snakemake.output.r2) + ".tmp"
target = int(snakemake.params.pairs)
written = 0

with open(snakemake.log[0], "w") as log, gzip.open(tmp1, "wt") as out1, gzip.open(tmp2, "wt") as out2:
    for r1_path, r2_path in zip(snakemake.input.r1, snakemake.input.r2):
        with gzip.open(r1_path, "rt") as in1, gzip.open(r2_path, "rt") as in2:
            while written < target:
                rec1 = [in1.readline() for _ in range(4)]
                rec2 = [in2.readline() for _ in range(4)]
                if not rec1[0] and not rec2[0]:
                    break
                if not all(rec1) or not all(rec2):
                    raise ValueError(f"incomplete paired FASTQ record in {r1_path} / {r2_path}")
                if normalize_read_id(rec1[0]) != normalize_read_id(rec2[0]):
                    raise ValueError(f"mate ID mismatch: {rec1[0].strip()} / {rec2[0].strip()}")
                out1.writelines(rec1)
                out2.writelines(rec2)
                written += 1
        if written >= target:
            break
    if written != target:
        raise ValueError(f"requested {target} fixture pairs but found {written}")
    log.write(f"sample={snakemake.wildcards.sample} fixture_pairs={written}\n")

os.replace(tmp1, snakemake.output.r1)
os.replace(tmp2, snakemake.output.r2)
atomic_json(snakemake.output.metadata, {
    "sample_id": str(snakemake.wildcards.sample),
    "pairs": written,
    "source_r1": list(snakemake.input.r1),
    "source_r2": list(snakemake.input.r2),
    "r1_sha256": sha256(snakemake.output.r1),
    "r2_sha256": sha256(snakemake.output.r2),
})
