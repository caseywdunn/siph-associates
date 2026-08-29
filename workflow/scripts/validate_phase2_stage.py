import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import ensure_parents

ensure_parents([snakemake.output[0], snakemake.log[0]])
errors = []
validations = [json.loads(Path(path).read_text()) for path in snakemake.input.validations]
expected = int(snakemake.params.expected_samples)

if len(validations) != expected:
    errors.append(f"validation count mismatch: {len(validations)} != {expected}")
if len({record.get("sample_id") for record in validations}) != len(validations):
    errors.append("duplicate sample validation IDs")
if any(record.get("status") != "PASS" for record in validations):
    errors.append("one or more sample validations did not pass")
if len({record.get("study") for record in validations}) != int(snakemake.params.expected_studies):
    errors.append("study coverage mismatch")
if len({record.get("host_route") for record in validations}) != int(snakemake.params.expected_routes):
    errors.append("host-route coverage mismatch")

with open(snakemake.input.qc, newline="") as handle:
    qc_rows = list(csv.DictReader(handle, delimiter="\t"))
if len(qc_rows) != expected:
    errors.append(f"QC aggregation count mismatch: {len(qc_rows)} != {expected}")
if {row.get("sample_id") for row in qc_rows} != {record.get("sample_id") for record in validations}:
    errors.append("QC aggregation sample IDs do not match validations")
if any(row.get("status") != "PASS" for row in qc_rows):
    errors.append("QC aggregation contains a non-PASS row")

with open(snakemake.input.nominations, newline="") as handle:
    nomination_reader = csv.DictReader(handle, delimiter="\t")
    required_nomination_columns = {"sample_id", "study", "source", "rank", "candidate_id", "candidate_name"}
    if not required_nomination_columns.issubset(set(nomination_reader.fieldnames or [])):
        errors.append("candidate nomination table lacks required columns")
    nomination_rows = list(nomination_reader)
if not nomination_rows:
    errors.append("candidate nomination aggregation contains no rows")
if not {row.get("sample_id") for row in nomination_rows}.issubset({record.get("sample_id") for record in validations}):
    errors.append("candidate nominations contain unexpected sample IDs")

provenance = json.loads(Path(snakemake.input.provenance).read_text())
if int(provenance.get("sample_count", -1)) != expected:
    errors.append("aggregation provenance sample count mismatch")

Path(snakemake.log[0]).write_text(
    "\n".join([f"samples={len(validations)}", f"nominations={len(nomination_rows)}", f"errors={len(errors)}", *errors]) + "\n"
)
if errors:
    raise ValueError("; ".join(errors))

temporary = Path(str(snakemake.output[0]) + ".tmp")
temporary.write_text(
    "status\tPASS\n"
    f"samples\t{len(validations)}\n"
    f"studies\t{len({record['study'] for record in validations})}\n"
    f"host_routes\t{len({record['host_route'] for record in validations})}\n"
    f"nominations\t{len(nomination_rows)}\n"
)
os.replace(temporary, snakemake.output[0])
