import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import ensure_parents

if not snakemake.params.enabled:
    raise ValueError("fixture-reference rule requested while fixture_enabled is false")

ensure_parents([snakemake.output[0], snakemake.log[0]])
temporary = str(snakemake.output[0]) + ".tmp"
limit = int(snakemake.params.bases)
written = 0
sequence_number = 0

with open(snakemake.input[0]) as source, open(temporary, "w") as target:
    sequence = []
    for line in source:
        if line.startswith(">"):
            if sequence and written < limit:
                chunk = "".join(sequence)[: limit - written]
                sequence_number += 1
                target.write(f">fixture_{snakemake.wildcards.route}_{sequence_number}\n")
                for offset in range(0, len(chunk), 80):
                    target.write(chunk[offset:offset + 80] + "\n")
                written += len(chunk)
            sequence = []
            if written >= limit:
                break
        else:
            sequence.append(line.strip())
    if sequence and written < limit:
        chunk = "".join(sequence)[: limit - written]
        sequence_number += 1
        target.write(f">fixture_{snakemake.wildcards.route}_{sequence_number}\n")
        for offset in range(0, len(chunk), 80):
            target.write(chunk[offset:offset + 80] + "\n")
        written += len(chunk)

if written != limit:
    raise ValueError(f"reference yielded {written} fixture bases, expected {limit}")
os.replace(temporary, snakemake.output[0])
Path(snakemake.log[0]).write_text(f"route={snakemake.wildcards.route} bases={written}\n")
