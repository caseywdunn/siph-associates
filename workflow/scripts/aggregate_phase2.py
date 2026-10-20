import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(
    0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts")
)
from common import atomic_json, ensure_parents, file_record, sha256
from screen_nominations import screen_nominations

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
sample_ids = list(snakemake.params.sample_ids)
validations = [
    json.loads(Path(path).read_text()) for path in snakemake.input.validations
]
records = {record["sample_id"]: record for record in validations}
if set(records) != set(sample_ids) or len(records) != len(sample_ids):
    raise ValueError("validation records do not exactly match aggregation sample IDs")
if any(record.get("status") != "PASS" for record in validations):
    raise ValueError("cannot aggregate non-PASS validation records")


def atomic_table(path, fields, rows):
    path = Path(str(path))
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


qc_fields = [
    "sample_id",
    "library_id",
    "specimen_id",
    "study",
    "host_route",
    "raw_pairs",
    "cap_pairs",
    "cap_applied",
    "pairs_entering_fastp",
    "pairs_after_fastp",
    "kraken_rows",
    "bracken_genus_rows",
    "bracken_species_rows",
    "sylph_rows",
    "phyloflash_outer_files",
    "status",
]
qc_rows = []
for sample in sample_ids:
    record = records[sample]
    counts = record["screen_rows"]
    qc_rows.append(
        {
            **record,
            "cap_applied": str(bool(record["cap_applied"])).lower(),
            "kraken_rows": counts["kraken"],
            "bracken_genus_rows": counts["genus"],
            "bracken_species_rows": counts["species"],
            "sylph_rows": counts["sylph"],
            "phyloflash_outer_files": counts["phyloflash_outer_files"],
        }
    )

nomination_fields, nominations = screen_nominations(
    sample_ids,
    {k: v["study"] for k, v in records.items()},
    snakemake.input.genus,
    snakemake.input.species,
    snakemake.input.sylph,
    snakemake.input.phyloflash,
)

atomic_table(snakemake.output.qc, qc_fields, qc_rows)
atomic_table(snakemake.output.nominations, nomination_fields, nominations)
atomic_json(
    snakemake.output.provenance,
    {
        "scope": str(snakemake.wildcards.scope),
        "sample_count": len(sample_ids),
        "sample_ids": sample_ids,
        "validation_sha256": {
            sample: sha256(path)
            for sample, path in zip(sample_ids, snakemake.input.validations)
        },
        "outputs": [
            file_record(snakemake.output.qc, checksum=True),
            file_record(snakemake.output.nominations, checksum=True),
        ],
    },
)
Path(snakemake.log[0]).write_text(
    f"scope={snakemake.wildcards.scope}\nsamples={len(sample_ids)}\nnominations={len(nominations)}\n"
)
