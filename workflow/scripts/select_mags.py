#!/usr/bin/env python3
"""Select MAGs meeting the locked CheckM2 quality rule and extract their FASTA files."""
import argparse
import csv
import json
import os
import tarfile
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--membership", required=True)
parser.add_argument("--results", required=True, type=Path)
parser.add_argument("--config", required=True)
parser.add_argument("--genomes", required=True, type=Path)
parser.add_argument("--table", required=True, type=Path)
args = parser.parse_args()

rules = json.loads(Path(args.config).read_text())["mags"]
fields = ["mag_id", "sample_id", "completeness", "contamination", "quality_score", "near_complete",
          "genome_size", "gc_content", "contig_n50", "checkm2_model"]
with open(args.membership, newline="") as handle:
    samples = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")
               if row["assembly_eligible"] == "true"]

args.genomes.mkdir(parents=True)
selected = []
for sample in samples:
    status = (args.results / "checkm2" / f"{sample}.status").read_text().strip()
    if status not in rules["checkm2_statuses_eligible"]:
        continue
    with open(args.results / "checkm2" / f"{sample}.tsv", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            completeness, contamination = float(row["Completeness"]), float(row["Contamination"])
            if completeness < rules["min_completeness"] or contamination >= rules["max_contamination_exclusive"]:
                continue
            selected.append({
                "mag_id": row["Name"], "sample_id": sample,
                "completeness": completeness, "contamination": contamination,
                "quality_score": round(completeness - 5 * contamination, 2),
                "near_complete": str(completeness >= rules["near_complete_min_completeness"]
                                     and contamination < rules["near_complete_max_contamination_exclusive"]).lower(),
                "genome_size": row["Genome_Size"], "gc_content": row["GC_Content"],
                "contig_n50": row["Contig_N50"], "checkm2_model": row["Completeness_Model_Used"],
            })
    wanted = {row["mag_id"] for row in selected if row["sample_id"] == sample}
    with tarfile.open(args.results / "bins" / f"{sample}.tar.gz") as archive:
        for member in archive.getmembers():
            if member.name.endswith(".fa") and member.name[:-3] in wanted:
                (args.genomes / member.name).write_bytes(archive.extractfile(member).read())
                wanted.discard(member.name[:-3])
    if wanted:
        raise SystemExit(f"CheckM2 bins missing from the archive for {sample}: {sorted(wanted)}")

temporary = args.table.with_name(args.table.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(selected)
os.replace(temporary, args.table)
print(f"selected {len(selected)} MAGs from {len({row['sample_id'] for row in selected})} libraries")
