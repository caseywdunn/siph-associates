import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import ensure_parents

ensure_parents([snakemake.output[0], snakemake.log[0]])
errors = []
with open(snakemake.input.samples, newline="") as handle:
    samples = list(csv.DictReader(handle, delimiter="\t"))
included = [row for row in samples if row["include_primary"].lower() == "true"]
if len(included) != int(snakemake.params.expected_samples):
    errors.append(f"included sample count {len(included)} != {snakemake.params.expected_samples}")
if len({row["sample_id"] for row in included}) != len(included):
    errors.append("sample_id values are not unique")
for row in included:
    r1 = [path for path in row["r1_paths"].split(";") if path]
    r2 = [path for path in row["r2_paths"].split(";") if path]
    if not r1 or len(r1) != len(r2):
        errors.append(f"{row['sample_id']}: missing or mismatched mate lists")
    for path in r1 + r2:
        if not Path(path).is_file():
            errors.append(f"{row['sample_id']}: missing {path}")
if not Path(snakemake.input.freeze).is_file():
    errors.append("manifest freeze is missing")

with open(snakemake.log[0], "w") as log:
    log.write(f"samples={len(included)} errors={len(errors)}\n")
    for error in errors:
        log.write("ERROR " + error + "\n")
if errors:
    raise ValueError("; ".join(errors[:10]))
temporary = str(snakemake.output[0]) + ".tmp"
Path(temporary).write_text(f"status\tPASS\nsamples\t{len(included)}\n")
os.replace(temporary, snakemake.output[0])
