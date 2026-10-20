#!/usr/bin/env python3
"""Export reviewed catalog tables and audit evidence, preserving original catalogs."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    out = Path("figures/catalog_review")
    out.mkdir(parents=True, exist_ok=True)
    exclusion_file = Path("config/catalog_domain_exclusions.tsv")
    exclusions = read(exclusion_file)
    targets = {r["catalog_id"]: r for r in exclusions}
    bins = {bid: r for r in exclusions for bid in r["bin_ids"].split(";")}
    sources = [exclusion_file, Path(__file__)]
    for name, source, key, rejected, field in [
        (
            "bacterial_catalog.manifest.tsv",
            "data/results/phase4_catalog/v2/bacterial_catalog.manifest.tsv",
            "catalog_id",
            targets,
            "status",
        ),
        (
            "mag_catalog.tsv",
            "data/results/phase3_catalog/mags/mag_catalog.tsv",
            "mag_id",
            bins,
            "catalog_status",
        ),
    ]:
        sources.append(Path(source))
        rows = read(source)
        assert set(rejected) <= {r[key] for r in rows}
        for row in rows:
            review = rejected.get(row[key])
            row["domain_review"] = "excluded_from_prokaryotes" if review else ""
            row["domain_review_reason"] = review["reason"] if review else ""
            row["reviewed_affinity"] = review["reviewed_affinity"] if review else ""
            if review:
                row[field] = "excluded_domain_review"
        write(out / name, rows)

    audit = Path("data/results/mag_domain_audit/2026-10-18/review")
    counts = read(audit / "genes_nr.counts.tsv")
    markers = read(audit / "genes_nr.archaeal_markers.tsv")
    genes = read(audit / "genes_nr.tsv")
    sources.extend(
        audit / n
        for n in [
            "genes_nr.counts.tsv",
            "genes_nr.archaeal_markers.tsv",
            "genes_nr.tsv",
        ]
    )
    rep = "Church2025__YPM-IZ-110469__bin.2"
    values = {
        "RepresentativeProteins": sum(
            int(r["queries"]) for r in counts if r["bin_id"] == rep
        ),
        "EukaryoticProteins": sum(
            int(r["queries"])
            for r in counts
            if r["bin_id"] == rep and r["assignment"] == "Eukaryota"
        ),
        "ArchaealProteins": sum(
            int(r["queries"])
            for r in counts
            if r["bin_id"] == rep and r["assignment"] == "Archaea"
        ),
        "PerkinselaProteins": sum(
            r["bin_id"] == rep
            and r["near_best_organisms"] == "Perkinsela sp. CCAP 1560/4"
            for r in genes
        ),
        "EukaryoticMarkers": sum(
            r["bin_id"] == rep and r["assignment"] == "Eukaryota" for r in markers
        ),
        "ExcludedBins": len(bins),
    }
    retained = [
        r for r in read(out / "mag_catalog.tsv") if r["catalog_status"] == "retained"
    ]
    values["RetainedMags"] = len(retained)
    values["MagClusters"] = len({r["species_cluster"] for r in retained})
    values["NearCompleteMags"] = sum(r["near_complete"] == "true" for r in retained)
    write(out / "gene_domain_counts.tsv", counts)
    write(
        out / "marker_domain_counts.tsv",
        [
            {"bin_id": b, "assignment": a, "hits": n}
            for (b, a), n in sorted(
                Counter((r["bin_id"], r["assignment"]) for r in markers).items()
            )
        ],
    )
    (out / "summary.json").write_text(json.dumps(values, indent=2) + "\n")
    (out / "values.tex").write_text(
        "".join(f"\\newcommand{{\\domain{k}}}{{{v:,}}}\n" for k, v in values.items())
    )
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    (out / "provenance.json").write_text(
        json.dumps(
            {
                "inputs_sha256": {str(p): sha(p) for p in sources},
                "outputs_sha256": {
                    p.name: sha(p) for p in out.iterdir() if p.name != "provenance.json"
                },
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(values, indent=2))


if __name__ == "__main__":
    main()
