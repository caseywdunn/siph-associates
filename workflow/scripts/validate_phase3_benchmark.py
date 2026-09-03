import csv
import json
import os
from pathlib import Path

table = list(csv.DictReader(open(snakemake.input.table), delimiter="\t"))
errors = []
if len(table) != int(snakemake.params.expected_rows):
    errors.append(f"expected {snakemake.params.expected_rows} rows, found {len(table)}")
if len({row["benchmark_id"] for row in table}) != len(table):
    errors.append("benchmark IDs are not unique")
if len({row["sample_id"] for row in table}) != int(snakemake.params.expected_samples):
    errors.append("sample count mismatch")
required = ("contigs", "assembly_bases", "rrna_markers", "viral_contigs", "raw_bins",
            "assessed_bins", "unassessed_bins")
for row in table:
    for field in required:
        try:
            if float(row[field]) < 0:
                errors.append(f"negative {field}: {row['benchmark_id']}")
        except Exception:
            errors.append(f"invalid {field}: {row['benchmark_id']}")
    if row.get("checkm2_status") not in {"assessed", "no_bins", "no_annotations"}:
        errors.append(f"invalid checkm2_status: {row['benchmark_id']}")
    try:
        if int(row["assessed_bins"]) + int(row["unassessed_bins"]) != int(row["raw_bins"]):
            errors.append(f"CheckM2 bin accounting mismatch: {row['benchmark_id']}")
    except Exception:
        pass
if {row["strategy"] for row in table} != {"reference_depleted", "kraken_nominated", "fixed_effort_trimmed"}:
    errors.append("strategy set mismatch")
provenance = json.loads(Path(snakemake.input.provenance).read_text())
if provenance.get("benchmark_rows") != len(table):
    errors.append("aggregation provenance row count mismatch")

output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
temporary = output.with_name(output.name + ".tmp")
lines = [
    f"status\t{'PASS' if not errors else 'FAIL'}",
    f"benchmark_rows\t{len(table)}",
    f"samples\t{len({row['sample_id'] for row in table})}",
    f"strategies\t{len({row['strategy'] for row in table})}",
    f"errors\t{len(errors)}",
]
lines.extend(f"ERROR\t{error}" for error in errors)
temporary.write_text("\n".join(lines) + "\n")
os.replace(temporary, output)
Path(snakemake.log[0]).write_text("\n".join(lines) + "\n")
if errors:
    raise ValueError("; ".join(errors))
