#!/usr/bin/env python3
"""Compare expanded publication outputs with the preserved pre-expansion baseline."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow/scripts"))
from eukaryote_units import NEGATIVE_CONTROLS, unit as reporting_unit  # noqa: E402

BASE = Path(
    "data/results/post_expansion/2026-10-21/reporting/baseline/figures-44ade17.tar.gz"
)
OUT = Path("figures/expanded_review")
FOCAL = {
    "MAGSP0005",
    "MAGSP0007",
    "MAGSP0010",
    "MAGSP0029",
    "MAGSP0031",
    "MAGSP0011",
    "MAGSP0012",
}


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name, rows, fields):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def background_flags(metadata, lineages, min_identity, min_pairs):
    """Flag libraries where a predefined background lineage reaches the read cutoff."""
    counts = defaultdict(Counter)
    for row in lineages:
        if int(row["min_identity_band"]) < min_identity:
            continue
        group = reporting_unit(row["lineage"])
        if group and group[0] == "negative_control":
            counts[row["sample_id"]][group[1]] += int(row["read_pairs"])
    results = []
    for row in metadata:
        sid = row["analysis_sample_id"]
        maximum = max(counts[sid].values(), default=0)
        results.append(
            {
                "sample_id": sid,
                "study": row["study"],
                "species_current": row["species_current"],
                **{f"{group}_pairs": counts[sid][group] for group in NEGATIVE_CONTROLS},
                "max_control_pairs": maximum,
                "min_identity_percent": min_identity,
                "primary_read_threshold": min_pairs,
                "background_flag": str(maximum >= min_pairs).lower(),
            }
        )
    return results


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with tarfile.open(BASE) as archive:

        def baseline(path):
            return list(
                csv.DictReader(
                    io.TextIOWrapper(archive.extractfile(path)), delimiter="\t"
                )
            )

        old_metadata = baseline("figures/cohort/sample_metadata.tsv")
        old_bacteria = baseline("figures/cohort/bacterial_detections.tsv")
        old_viruses = baseline("figures/cohort/viral_detections.tsv")
        old_eukaryotes = baseline("figures/cohort/eukaryote_evidence.tsv")
        old_catalog = baseline("figures/catalog_review/bacterial_catalog.manifest.tsv")
        old_mags = baseline("figures/catalog_review/mag_catalog.tsv")
        old_summary = json.load(archive.extractfile("figures/cohort/summary.json"))
    new_metadata = read("figures/cohort/sample_metadata.tsv")
    settings = json.loads(Path("config/reporting_sources.json").read_text())
    presence = json.loads(Path("config/eukaryote_gate.json").read_text())["presence"][
        "validated_reads"
    ]
    flags = background_flags(
        new_metadata,
        read(Path(settings["eukaryote_root"]) / "read_lineages.tsv"),
        presence["min_identity_band"],
        presence["min_read_pairs"],
    )
    preservation = {
        r["analysis_sample_id"]: r["preservation"]
        for r in read("data/metadata/manko2026/sample_metadata.tsv")
    }
    for row in flags:
        row["recorded_preservation"] = preservation.get(row["sample_id"], "")
    write("ssu_background_flags.tsv", flags, list(flags[0]))
    flagged = {r["sample_id"] for r in flags if r["background_flag"] == "true"}
    original = {r["analysis_sample_id"] for r in old_metadata}
    expected = {r["sample_id"] for r in read("config/post_expansion_samples.tsv")}
    if {r["analysis_sample_id"] for r in new_metadata} != expected or len(
        new_metadata
    ) != len(expected):
        raise ValueError("Publication metadata do not match the expanded cohort")
    old_by_genome = {r["genome_id"]: r for r in old_catalog}
    catalog = read("figures/catalog_qc/catalog_v3/bacterial_catalog.manifest.tsv")
    identical = []
    for row in catalog:
        if row["genome_id"] not in FOCAL:
            continue
        previous = old_by_genome[row["genome_id"]]
        same = all(
            row[k] == previous[k] for k in ["origin", "origin_sha256", "fasta_sha256"]
        )
        identical.append(
            {
                "genome_id": row["genome_id"],
                "catalog_id": row["catalog_id"],
                "identical_source_and_sequence": str(same).lower(),
                "origin_sha256": row["origin_sha256"],
                "fasta_sha256": row["fasta_sha256"],
            }
        )
    if len(identical) != len(FOCAL) or any(
        r["identical_source_and_sequence"] != "true" for r in identical
    ):
        raise ValueError(
            "Existing focal phylogeny/annotations cannot be retained without review"
        )
    write("retained_focal_genomes.tsv", identical, list(identical[0]))
    metadata_changes = []
    previous_metadata = {r["analysis_sample_id"]: r for r in old_metadata}
    for row in new_metadata:
        sid = row["analysis_sample_id"]
        if sid not in previous_metadata:
            continue
        for field, value in previous_metadata[sid].items():
            if row.get(field, "") != value:
                metadata_changes.append(
                    dict(
                        analysis_sample_id=sid,
                        field=field,
                        previous=value,
                        current=row.get(field, ""),
                    )
                )
    write(
        "original_metadata_changes.tsv",
        metadata_changes,
        ["analysis_sample_id", "field", "previous", "current"],
    )
    crosswalk = read("figures/catalog_qc/catalog_v3/catalog_id_crosswalk.tsv")
    changes, totals, background_sensitivity = [], {}, []
    for domain, previous, current_path, unit in [
        (
            "microbial",
            old_bacteria,
            "figures/cohort/bacterial_detections.tsv",
            "target_id",
        ),
        ("viral", old_viruses, "figures/cohort/viral_detections.tsv", "target_id"),
        (
            "eukaryotic",
            old_eukaryotes,
            "figures/cohort/eukaryote_evidence.tsv",
            "reporting_unit",
        ),
    ]:
        current = read(current_path)
        if domain == "eukaryotic":
            previous = [r for r in previous if r["count_as_detection"] == "true"]
            current = [r for r in current if r["count_as_detection"] == "true"]
        kind = (
            "bacterial"
            if domain == "microbial"
            else "viral"
            if domain == "viral"
            else ""
        )
        translation = {
            r["previous_id"]: r["current_id"]
            for r in crosswalk
            if r["kind"] == kind and r["previous_id"] and r["current_id"]
        }
        before, previous_ids = {}, defaultdict(set)
        if len({(r["sample_id"], r[unit]) for r in previous}) != len(previous):
            raise ValueError(f"Duplicate original {domain} detection units")
        for row in previous:
            key = row["sample_id"], translation.get(row[unit], row[unit])
            previous_ids[key].add(row[unit])
            if key not in before or row["grade"] == "high_confidence":
                before[key] = row["grade"]
        after = {(r["sample_id"], r[unit]): r["grade"] for r in current}
        if len(after) != len(current):
            raise ValueError(f"Duplicate {domain} detection units")
        for sid, target in sorted(set(before) | set(after)):
            a, b = before.get((sid, target), ""), after.get((sid, target), "")
            if a != b:
                change = (
                    "new_library"
                    if sid not in original
                    else "lost_detection"
                    if not b
                    else "gained_detection"
                    if not a
                    else "grade_change"
                )
                changes.append(
                    dict(
                        domain=domain,
                        sample_id=sid,
                        target=target,
                        previous_targets=";".join(sorted(previous_ids[sid, target])),
                        previous_grade=a,
                        current_grade=b,
                        change=change,
                    )
                )
        totals[domain] = {
            "previous_detections": len(previous),
            "previous_units_after_catalog_crosswalk": len(before),
            "current_detections": len(current),
            "current_in_original_libraries": sum(
                r["sample_id"] in original for r in current
            ),
            "current_in_added_libraries": sum(
                r["sample_id"] not in original for r in current
            ),
            "changes": dict(
                Counter(r["change"] for r in changes if r["domain"] == domain)
            ),
        }
        for scope, omitted in [
            ("all_libraries", set()),
            ("excluding_background_flagged_libraries", flagged),
        ]:
            retained = [r for r in current if r["sample_id"] not in omitted]
            background_sensitivity.append(
                {
                    "domain": domain,
                    "scope": scope,
                    "libraries_in_scope": len(expected - omitted),
                    "detections": len(retained),
                    "positive_libraries": len({r["sample_id"] for r in retained}),
                    "detected_targets_or_groups": len({r[unit] for r in retained}),
                    "high_confidence_detections": sum(
                        r["grade"] == "high_confidence" for r in retained
                    ),
                }
            )
    write(
        "detection_changes.tsv",
        changes,
        [
            "domain",
            "sample_id",
            "target",
            "previous_targets",
            "previous_grade",
            "current_grade",
            "change",
        ],
    )
    write(
        "background_sensitivity.tsv",
        background_sensitivity,
        list(background_sensitivity[0]),
    )
    old_retained = {r["mag_id"] for r in old_mags if r["catalog_status"] == "retained"}
    new_retained = {
        r["mag_id"]
        for r in read("figures/catalog_qc/catalog_v3/mag_catalog.tsv")
        if r["catalog_status"] == "retained"
    }
    summary = {
        "baseline_revision": "44ade17",
        "baseline_sha256": hashlib.sha256(BASE.read_bytes()).hexdigest(),
        "libraries": len(expected),
        "original_libraries": len(original),
        "added_libraries": len(expected - original),
        "original_metadata_changes": len(metadata_changes),
        "detection_comparison": totals,
        "ssu_background_flagged_libraries": sorted(flagged),
        "background_sensitivity": background_sensitivity,
        "background_interpretation": "Primary grades and thresholds are unchanged. The predefined SSU background-control lineages flag specimen-specific uncertainty. Omitting flagged libraries is a descriptive sensitivity check, not proof that their microbial or viral detections are false.",
        "retained_mag_ids_added": sorted(new_retained - old_retained),
        "retained_mag_ids_removed": sorted(old_retained - new_retained),
        "prior_cohort_summary": old_summary,
        "current_cohort_summary": json.loads(
            Path("figures/cohort/summary.json").read_text()
        ),
        "interpretation": "Catalog target IDs use the frozen v3 crosswalk. Changes include reference competition, catalog membership, assembly support and added libraries; read counts are not independent organisms.",
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in summary.items()
                if k not in ["prior_cohort_summary", "current_cohort_summary"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
