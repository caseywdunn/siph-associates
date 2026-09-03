import csv
import json
import os
from pathlib import Path

rows = list(csv.DictReader(open(snakemake.input.table), delimiter="\t"))
errors = []
if len(rows) != int(snakemake.params.expected_samples):
    errors.append(f"expected {snakemake.params.expected_samples} rows, found {len(rows)}")
if len({row["sample_id"] for row in rows}) != len(rows):
    errors.append("sample IDs are not unique")
if len({row["host_route"] for row in rows}) != int(snakemake.params.expected_routes):
    errors.append("host-route count mismatch")
expected = {"P_physalis": "reference_depleted", "N_septata": "reference_depleted",
            "none": "fixed_effort_trimmed"}
for row in rows:
    if expected.get(row["host_route"]) != row["strategy"]:
        errors.append(f"strategy mismatch: {row['sample_id']}")
    try:
        retained = int(row["retained_pairs"])
        threshold = int(row["minimum_assembly_input_pairs"])
        if (row["assembly_eligible"] == "true") != (retained >= threshold):
            errors.append(f"eligibility mismatch: {row['sample_id']}")
    except Exception:
        errors.append(f"invalid pair counts: {row['sample_id']}")
provenance = json.loads(Path(snakemake.input.provenance).read_text())
if provenance.get("samples") != len(rows):
    errors.append("provenance sample count mismatch")

output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
temporary = output.with_name(output.name + ".tmp")
lines = [f"status\t{'PASS' if not errors else 'FAIL'}", f"samples\t{len(rows)}",
         f"eligible\t{sum(row.get('assembly_eligible') == 'true' for row in rows)}",
         f"errors\t{len(errors)}", *[f"ERROR\t{error}" for error in errors]]
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, output)
Path(snakemake.log[0]).write_text("\n".join(lines) + "\n")
if errors:
    raise ValueError("; ".join(errors))
