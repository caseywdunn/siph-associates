import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import ensure_parents

ensure_parents([snakemake.output[0], snakemake.log[0]])
records = [json.loads(Path(path).read_text()) for path in snakemake.input]
errors = []
if len(records) != int(snakemake.params.expected_samples):
    errors.append("sample count mismatch")
if len({record["sample_id"] for record in records}) != len(records):
    errors.append("duplicate sample validations")
if any(record.get("status") != "PASS" for record in records):
    errors.append("one or more sample validations did not pass")
if len({record["study"] for record in records}) != int(snakemake.params.expected_studies):
    errors.append("study coverage mismatch")
if len({record["host_route"] for record in records}) != int(snakemake.params.expected_routes):
    errors.append("host-route coverage mismatch")
Path(snakemake.log[0]).write_text(f"samples={len(records)} errors={len(errors)}\n" + "\n".join(errors))
if errors:
    raise ValueError("; ".join(errors))
temporary = str(snakemake.output[0]) + ".tmp"
Path(temporary).write_text(
    "status\tPASS\n"
    f"samples\t{len(records)}\n"
    f"studies\t{len({r['study'] for r in records})}\n"
    f"host_routes\t{len({r['host_route'] for r in records})}\n"
)
os.replace(temporary, snakemake.output[0])
