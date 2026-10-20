#!/usr/bin/env python3
"""Combine domain, contig, coverage, coding-density and marker evidence for review."""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from assembly_completion import table, unique, write_table

LABELS = (
    "Bacteria",
    "Archaea",
    "Eukaryota",
    "Viruses",
    "ambiguous",
    "unresolved_taxonomy",
    "no_qualifying_hit",
)


def covered(intervals: list[tuple[int, int]]) -> int:
    """Count covered bases using half-open intervals without double counting overlaps."""
    total, end = 0, 0
    for start, stop in sorted(intervals):
        total += max(0, stop - max(start, end))
        end = max(end, stop)
    return total


def main() -> None:
    config = json.loads(Path("config/catalog_qc.json").read_text())
    root = Path(config["run"])
    prep = Path("data/results/post_expansion/2026-10-20")
    candidates = unique(table("figures/catalog_qc/candidate_review.tsv"), "mag_id")
    contigs = unique(table(prep / "qc/contig_review.tsv"), "contig_id")
    counts, by_contig, intervals = {}, {}, defaultdict(list)
    for source in ("genes", "tiles"):
        counts[source], by_contig[source] = defaultdict(Counter), defaultdict(Counter)
        with (root / f"review/{source}_swissprot.tsv").open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name, contig = row["bin_id"], row["contig_id"]
                if name not in candidates or contigs[contig]["mag_id"] != name:
                    raise ValueError(
                        f"Query/bin/contig identity mismatch: {row['query_id']}"
                    )
                start, end = int(row["start"]), int(row["end"])
                if not 1 <= start <= end <= int(contigs[contig]["length"]):
                    raise ValueError(
                        f"Gene/tile coordinate out of bounds: {row['query_id']}"
                    )
                counts[source][name][row["assignment"]] += 1
                by_contig[source][contig][row["assignment"]] += 1
                if source == "genes":
                    intervals[contig].append((start - 1, end))
    depth = {}
    for sample in sorted({r["sample_id"] for r in candidates.values()}):
        for row in table(prep / "phase3_cohort/bins" / f"{sample}.depth.tsv"):
            if row["contigName"] in contigs:
                depth[row["contigName"]] = float(row["totalAvgDepth"])
    markers = {}
    for kind in ("bac120", "ar53"):
        markers[kind] = unique(
            table(
                root / "gtdb_identify/identify" / f"gtdbtk.{kind}.markers_summary.tsv"
            ),
            "name",
        )
    out = Path("figures/catalog_qc")
    detail, bin_records = [], []
    for name, candidate in candidates.items():
        records = [r for r in contigs.values() if r["mag_id"] == name]
        coding = sum(covered(intervals[r["contig_id"]]) for r in records)
        bases = sum(int(r["length"]) for r in records)
        if bases != int(candidate["genome_size"]):
            raise ValueError(f"CheckM2 sequence size mismatch: {name}")
        row = dict(candidate, coding_fraction=coding / bases)
        for source in counts:
            row.update(
                {f"{source}_{label}": counts[source][name][label] for label in LABELS}
            )
        for kind in markers:
            for field in (
                "number_unique_genes",
                "number_multiple_genes",
                "number_multiple_unique_genes",
                "number_missing_genes",
            ):
                row[f"{kind}_{field}"] = markers[kind][name][field]
        row["mean_assembly_coverage"] = (
            sum(int(r["length"]) * depth[r["contig_id"]] for r in records) / bases
        )
        bin_records.append(row)
        for record in records:
            contig = record["contig_id"]
            item = dict(
                record,
                assembly_coverage=depth[contig],
                coding_fraction=covered(intervals[contig]) / int(record["length"]),
            )
            for source in counts:
                item.update(
                    {
                        f"{source}_{label}": by_contig[source][contig][label]
                        for label in LABELS
                    }
                )
            detail.append(item)
    write_table(out / "domain_review.tsv", bin_records, list(bin_records[0]))
    write_table(prep / "qc/contig_evidence.tsv", detail, list(detail[0]))
    print(json.dumps({"bins": len(bin_records), "contigs": len(detail)}, indent=2))


if __name__ == "__main__":
    main()
