import csv
import io
import json
import os
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, ensure_parents, file_record, sha256

ensure_parents(list(snakemake.output) + [snakemake.log[0]])
sample_ids = list(snakemake.params.sample_ids)
validations = [json.loads(Path(path).read_text()) for path in snakemake.input.validations]
records = {record["sample_id"]: record for record in validations}
if set(records) != set(sample_ids) or len(records) != len(sample_ids):
    raise ValueError("validation records do not exactly match aggregation sample IDs")
if any(record.get("status") != "PASS" for record in validations):
    raise ValueError("cannot aggregate non-PASS validation records")


def atomic_table(path, fields, rows):
    path = Path(str(path))
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


qc_fields = [
    "sample_id", "library_id", "specimen_id", "study", "host_route", "raw_pairs",
    "cap_pairs", "cap_applied", "pairs_entering_fastp", "pairs_after_fastp",
    "kraken_rows", "bracken_genus_rows", "bracken_species_rows", "sylph_rows",
    "phyloflash_outer_files", "status",
]
qc_rows = []
for sample in sample_ids:
    record = records[sample]
    counts = record["screen_rows"]
    qc_rows.append({
        **record,
        "cap_applied": str(bool(record["cap_applied"])).lower(),
        "kraken_rows": counts["kraken"],
        "bracken_genus_rows": counts["genus"],
        "bracken_species_rows": counts["species"],
        "sylph_rows": counts["sylph"],
        "phyloflash_outer_files": counts["phyloflash_outer_files"],
    })

nomination_fields = [
    "sample_id", "study", "source", "rank", "candidate_id", "candidate_name",
    "estimated_reads", "fraction_total_reads", "taxonomic_abundance",
    "sequence_abundance", "adjusted_ani", "true_coverage", "raw_support",
]
nominations = []

for sample, genus_path, species_path, sylph_path, phylo_path in zip(
        sample_ids, snakemake.input.genus, snakemake.input.species,
        snakemake.input.sylph, snakemake.input.phyloflash):
    study = records[sample]["study"]
    for rank, path in (("genus", genus_path), ("species", species_path)):
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                nominations.append({
                    "sample_id": sample, "study": study, "source": "Bracken", "rank": rank,
                    "candidate_id": row["taxonomy_id"], "candidate_name": row["name"],
                    "estimated_reads": row["new_est_reads"],
                    "fraction_total_reads": row["fraction_total_reads"],
                    "raw_support": row["kraken_assigned_reads"],
                })
    with open(sylph_path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            nominations.append({
                "sample_id": sample, "study": study, "source": "sylph", "rank": "genome",
                "candidate_id": row["Genome_file"], "candidate_name": Path(row["Genome_file"]).name,
                "taxonomic_abundance": row["Taxonomic_abundance"],
                "sequence_abundance": row["Sequence_abundance"],
                "adjusted_ani": row["Adjusted_ANI"], "true_coverage": row["True_cov"],
                "raw_support": row["kmers_reassigned"],
            })
    with tarfile.open(phylo_path, "r:gz") as outer:
        nested_member = next(member for member in outer.getmembers() if member.name.endswith(".phyloFlash.tar.gz"))
        nested_stream = outer.extractfile(nested_member)
        if nested_stream is None:
            raise ValueError(f"cannot read nested phyloFlash archive for {sample}")
        with tarfile.open(fileobj=io.BytesIO(nested_stream.read()), mode="r:gz") as inner:
            ntu_member = next(member for member in inner.getmembers() if member.name.endswith(".phyloFlash.NTUfull_abundance.csv"))
            ntu_stream = inner.extractfile(ntu_member)
            if ntu_stream is None:
                raise ValueError(f"cannot read phyloFlash NTU table for {sample}")
            for line in io.TextIOWrapper(ntu_stream, encoding="utf-8"):
                lineage, count = line.rstrip("\n").rsplit(",", 1)
                nominations.append({
                    "sample_id": sample, "study": study, "source": "phyloFlash",
                    "rank": "SSU_lineage", "candidate_id": lineage,
                    "candidate_name": lineage.split(";")[-1], "raw_support": count,
                })

atomic_table(snakemake.output.qc, qc_fields, qc_rows)
atomic_table(snakemake.output.nominations, nomination_fields, nominations)
atomic_json(snakemake.output.provenance, {
    "scope": str(snakemake.wildcards.scope),
    "sample_count": len(sample_ids),
    "sample_ids": sample_ids,
    "validation_sha256": {sample: sha256(path) for sample, path in zip(sample_ids, snakemake.input.validations)},
    "outputs": [file_record(snakemake.output.qc, checksum=True), file_record(snakemake.output.nominations, checksum=True)],
})
Path(snakemake.log[0]).write_text(
    f"scope={snakemake.wildcards.scope}\nsamples={len(sample_ids)}\nnominations={len(nominations)}\n"
)
