import csv
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(snakemake.config["manifest"]).resolve().parent / "workflow" / "scripts"))
from common import atomic_json, file_record, sha256


def records(paths):
    return {record["benchmark_id"]: record for record in (json.loads(Path(path).read_text()) for path in paths)}


def benchmark_record(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty Snakemake benchmark: {path}")
    return rows[-1]


panel = list(csv.DictReader(open(snakemake.input.panel), delimiter="\t"))
ids = [row["benchmark_id"] for row in panel]
assemblies = records(snakemake.input.assemblies)
markers = records(snakemake.input.markers)
viruses = records(snakemake.input.viruses)
bins = records(snakemake.input.bins)
checkm2 = records(snakemake.input.checkm2)
host = records(snakemake.input.host)
for name, mapping in (("assemblies", assemblies), ("markers", markers), ("viruses", viruses),
                      ("bins", bins), ("checkm2", checkm2), ("host", host)):
    if set(mapping) != set(ids):
        raise ValueError(f"{name} metrics do not match benchmark panel")

resource_sets = {}
for label, paths in (("assembly", snakemake.input.assembly_bench), ("virus", snakemake.input.virus_bench),
                     ("binning", snakemake.input.bin_bench), ("checkm2", snakemake.input.checkm2_bench)):
    resource_sets[label] = {Path(path).stem: benchmark_record(path) for path in paths}

fields = [
    "benchmark_id", "sample_id", "strategy", "signal_stratum", "selection_basis",
    "contigs", "assembly_bases", "n50", "max_contig", "rrna_markers", "bacterial_rrna",
    "archaeal_rrna", "eukaryotic_rrna", "viral_contigs", "viral_bases", "raw_bins",
    "binned_bases", "assessed_bins", "medium_quality_or_better_bins", "high_quality_bins",
    "host_mapped_record_fraction", "assembly_seconds", "assembly_max_rss_mb", "virus_seconds",
    "virus_max_rss_mb", "binning_seconds", "binning_max_rss_mb", "checkm2_seconds",
    "checkm2_max_rss_mb",
]
rows = []
for panel_row in panel:
    key = panel_row["benchmark_id"]
    def resource(label, field):
        return resource_sets[label][key].get(field, "")
    rows.append({
        **panel_row,
        "contigs": assemblies[key]["contigs"], "assembly_bases": assemblies[key]["total_bases"],
        "n50": assemblies[key]["n50"], "max_contig": assemblies[key]["max_contig"],
        "rrna_markers": markers[key]["total_rrna"], "bacterial_rrna": markers[key]["bacterial_rrna"],
        "archaeal_rrna": markers[key]["archaeal_rrna"], "eukaryotic_rrna": markers[key]["eukaryotic_rrna"],
        "viral_contigs": viruses[key]["viral_contigs"], "viral_bases": viruses[key]["viral_bases"],
        "raw_bins": bins[key]["raw_bins"], "binned_bases": bins[key]["binned_bases"],
        "assessed_bins": checkm2[key]["assessed_bins"],
        "medium_quality_or_better_bins": checkm2[key]["medium_quality_or_better_bins"],
        "high_quality_bins": checkm2[key]["high_quality_bins"],
        "host_mapped_record_fraction": host[key].get("host_mapped_record_fraction"),
        "assembly_seconds": resource("assembly", "s"), "assembly_max_rss_mb": resource("assembly", "max_rss"),
        "virus_seconds": resource("virus", "s"), "virus_max_rss_mb": resource("virus", "max_rss"),
        "binning_seconds": resource("binning", "s"), "binning_max_rss_mb": resource("binning", "max_rss"),
        "checkm2_seconds": resource("checkm2", "s"), "checkm2_max_rss_mb": resource("checkm2", "max_rss"),
    })

output = Path(snakemake.output.table)
output.parent.mkdir(parents=True, exist_ok=True)
temporary = output.with_name(output.name + ".tmp")
with temporary.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
os.replace(temporary, output)
atomic_json(snakemake.output.provenance, {
    "benchmark_rows": len(rows),
    "samples": sorted({row["sample_id"] for row in rows}),
    "strategies": sorted({row["strategy"] for row in rows}),
    "panel_sha256": sha256(snakemake.input.panel),
    "config_sha256": sha256(snakemake.input.config),
    "output": file_record(output, checksum=True),
})
Path(snakemake.log[0]).write_text(f"rows={len(rows)} samples={len(set(row['sample_id'] for row in rows))}\n")
