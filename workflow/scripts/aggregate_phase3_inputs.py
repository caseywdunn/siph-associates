import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, file_record, sha256

sample_ids = list(snakemake.params.sample_ids)
records = [json.loads(Path(path).read_text()) for path in snakemake.input.metrics]
by_sample = {record["sample_id"]: record for record in records}
if len(records) != len(by_sample) or set(by_sample) != set(sample_ids):
    raise ValueError("Phase-3 input metrics do not exactly match the cohort")

fields = ["sample_id", "host_route", "strategy", "input_pairs", "retained_pairs",
          "retained_fraction", "minimum_assembly_input_pairs", "assembly_eligible"]
rows = []
for sample in sample_ids:
    record = by_sample[sample]
    eligible = int(record["retained_pairs"]) >= int(snakemake.params.minimum_pairs)
    rows.append({**record, "minimum_assembly_input_pairs": int(snakemake.params.minimum_pairs),
                 "assembly_eligible": str(eligible).lower()})

output = Path(snakemake.output.table)
output.parent.mkdir(parents=True, exist_ok=True)
temporary = output.with_name(output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
os.replace(temporary, output)
atomic_json(snakemake.output.provenance, {
    "samples": len(rows),
    "eligible_samples": sum(row["assembly_eligible"] == "true" for row in rows),
    "minimum_assembly_input_pairs": int(snakemake.params.minimum_pairs),
    "config_sha256": sha256(snakemake.input.config),
    "output": file_record(output, checksum=True),
})
Path(snakemake.log[0]).write_text(
    f"samples={len(rows)}\neligible={sum(row['assembly_eligible'] == 'true' for row in rows)}\n"
)
