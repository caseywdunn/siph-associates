#!/usr/bin/env python3
"""Describe accepted associate detections by collection depth, host and tissue.

All libraries contribute to descriptive summaries. Single-flowcell subsets are
additional descriptive checks; this script fits no models and computes no tests.
Depth intervals remain intervals, and missing detections mean no accepted signal
in the supplied evidence tables rather than demonstrated biological absence.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import shlex
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

VALID = {"validated", "high_confidence"}
COUNTS = {
    "bacterial": "bacterial_targets",
    "archaeal": "archaeal_targets",
    "viral": "viral_targets",
    "eukaryote": "eukaryote_groups",
    "dt68": "dt68_targets",
    "vibrio": "vibrio_targets",
    "metamycoplasmataceae": "metamycoplasmataceae_targets",
    "clade_a": "clade_a_targets",
    "clade_b": "clade_b_targets",
    "prey": "prey_groups",
    "parasite": "parasite_groups",
}
FOCAL = ("MAGSP0005", "MAGSP0007", "MAGSP0010", "Vibrio")
INPUTS = {
    "manifest": "manifest.csv",
    "bacteria": "figures/cohort/bacterial_detections.tsv",
    "viruses": "figures/cohort/viral_detections.tsv",
    "eukaryotes": "figures/cohort/eukaryote_evidence.tsv",
    "flowcells": "figures/cohort/library_flowcells.tsv",
    "lineages": "figures/cohort/mycoplasmatales_incidence.tsv",
}


def read_table(path: Path) -> list[dict[str, str]]:
    """Read small reviewed CSV/TSV inputs without changing values."""
    with path.open(newline="") as handle:
        return list(
            csv.DictReader(handle, delimiter="," if path.suffix == ".csv" else "\t")
        )


def write_table(path: Path, rows: list[dict]) -> None:
    """Write deterministic human-readable tables."""
    if not rows:
        raise ValueError(f"No rows for required output {path.name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def index_manifest(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Require a bijection between canonical and analysis sample identifiers."""
    output = {}
    canonical = set()
    for row in rows:
        library = row["library_id"]
        if library.count(":") != 1 or "__" in library or library in canonical:
            raise ValueError(f"Ambiguous or duplicate library ID: {library}")
        sample = library.replace(":", "__", 1)
        if sample in output:
            raise ValueError(f"Duplicate analysis ID: {sample}")
        canonical.add(library)
        output[sample] = row
    return output


def unique_evidence(
    rows: list[dict[str, str]], target: str, meta: dict
) -> list[dict[str, str]]:
    """Collapse exact duplicate evidence only; conflicting copies are errors."""
    seen = {}
    for row in rows:
        sample = row["sample_id"]
        if sample not in meta:
            raise ValueError(f"Unmatched evidence sample: {sample}")
        key = (sample, row[target])
        if key in seen and seen[key] != row:
            raise ValueError(f"Conflicting evidence rows: {key}")
        seen[key] = row
    return [seen[key] for key in sorted(seen)]


def numeric(value: str) -> float | None:
    """Validate nonnegative finite collection depths, preserving missingness."""
    if not value:
        return None
    number = float(value)
    if not math.isfinite(number) or number < 0:
        raise ValueError(f"Invalid depth: {value}")
    return number


def depth_context(row: dict[str, str]) -> dict[str, str]:
    """Retain recorded point/range and classify only unambiguous depth strata."""
    point = numeric(row.get("depth_m", ""))
    lower = numeric(row.get("depth_min_m", ""))
    upper = numeric(row.get("depth_max_m", ""))
    if (lower is None) != (upper is None):
        raise ValueError(f"Incomplete depth interval: {row['library_id']}")
    if lower is not None and (
        lower > upper or (point is not None and not lower <= point <= upper)
    ):
        raise ValueError(f"Inconsistent depth interval: {row['library_id']}")
    assigned = row.get("collection_depth_basis") == "curator_assigned_surface"
    if assigned and (point != 0 or not row["species_current"].startswith("Physalia ")):
        raise ValueError("Surface assignment must be a zero-depth Physalia record")
    if point is None and lower is None:
        status = "unknown"
    elif assigned:
        status = "assigned_surface"
    elif lower is not None:
        status = "recorded_point_and_range" if point is not None else "recorded_range"
    else:
        status = "recorded_point"
    # The full interval, when supplied, determines the stratum; no midpoint is imputed.
    lo, hi = (lower, upper) if lower is not None else (point, point)
    if row["species_current"].startswith("Physalia "):
        if lo != 0 or hi != 0:
            raise ValueError("Expected accepted surface convention for all Physalia")
        group = "Physalia_surface"
    elif lo is None:
        group = "non_Physalia_unknown"
    elif hi <= 20:
        group = "non_Physalia_0_20m"
    elif lo > 20:
        group = "non_Physalia_gt20m"
    else:
        group = "non_Physalia_crosses20m"
    return {"depth_status": status, "depth_group": group}


def lineage_mapping(rows: list[dict[str, str]]) -> dict[str, tuple[str, str]]:
    """Read the already reviewed catalog-to-lineage mapping from cohort export."""
    mapping = {}
    genomes = {}
    for row in rows:
        target = row["catalog_id"]
        pair = (row["species_cluster"], row["lineage"])
        if target in mapping and mapping[target] != pair:
            raise ValueError(f"Conflicting focal lineage assignment: {target}")
        if pair[0] in genomes and genomes[pair[0]] != target:
            raise ValueError(f"Nonunique focal genome assignment: {pair[0]}")
        mapping[target] = pair
        genomes[pair[0]] = target
    required = set(FOCAL[:-1]) | {"MAGSP0011"}
    if not required.issubset(genomes):
        raise ValueError("Required focal targets missing from reviewed lineage table")
    return mapping


def derive_context(
    inputs: dict[str, list[dict[str, str]]],
) -> tuple[list[dict], list[dict]]:
    """Join all accepted evidence to every manifest library and retain raw support."""
    meta = index_manifest(inputs["manifest"])
    mapping = lineage_mapping(inputs["lineages"])
    flow_rows = unique_evidence(inputs["flowcells"], "sample_id", meta)
    flow = {row["sample_id"]: row for row in flow_rows}
    if set(flow) != set(meta):
        raise ValueError("Flowcell table must cover exactly the manifest libraries")
    present = {sample: defaultdict(set) for sample in meta}
    grades = defaultdict(list)
    evidence = []
    for domain, target_column, key in [
        ("prokaryote", "target_id", "bacteria"),
        ("viral", "target_id", "viruses"),
        ("eukaryote", "reporting_unit", "eukaryotes"),
    ]:
        rows = unique_evidence(inputs[key], target_column, meta)
        for row in rows:
            selected = row["grade"] in VALID
            if key == "bacteria":
                selected &= row["role"] != "decoy"
            elif key == "viruses":
                selected &= row["catalog_class"] == "associate"
            else:
                if row["count_as_detection"] not in {"true", "false"}:
                    raise ValueError("Eukaryote reporting flag must be true/false")
                if row["count_as_detection"] == "true" and not selected:
                    raise ValueError("Eukaryote accepted flag conflicts with grade")
                selected &= row["count_as_detection"] == "true"
            if not selected:
                continue
            sample, target = row["sample_id"], row[target_column]
            tax = set(row.get("taxonomy", "").split(";"))
            actual_domain = domain
            if domain == "prokaryote":
                if "d__Archaea" in tax:
                    actual_domain = "archaeal"
                elif "d__Bacteria" in tax:
                    actual_domain = "bacterial"
                else:
                    raise ValueError(
                        f"Unrecognized prokaryote domain: {row['taxonomy']}"
                    )
            present[sample][actual_domain].add(target)
            genome, lineage = mapping.get(target, ("", ""))
            if "g__Vibrio" in tax:
                present[sample]["vibrio"].add(target)
            if "f__Metamycoplasmataceae" in tax:
                present[sample]["metamycoplasmataceae"].add(target)
            for label, bucket in (
                ("Clade A", "clade_a"),
                ("Clade B", "clade_b"),
                ("DT-68", "dt68"),
            ):
                if lineage == label:
                    present[sample][bucket].add(target)
            if lineage == "DT-68":
                grades[sample].append(row["grade"])
            if genome in FOCAL:
                present[sample][genome].add(target)
            role = row.get("role", "") if domain == "eukaryote" else ""
            if role in {"prey", "parasite"}:
                present[sample][role].add(target)
            evidence.append(
                {
                    "library_id": meta[sample]["library_id"],
                    "sample_id": sample,
                    "domain": actual_domain,
                    "target_or_group": target,
                    "unit": "eukaryote_reporting_group"
                    if domain == "eukaryote"
                    else "genome_target",
                    "grade": row["grade"],
                    "genome_id": genome,
                    "focal_lineage": lineage,
                    "taxonomy": row.get("taxonomy", row.get("named_lineage", "")),
                    "biological_role": role,
                    "breadth": row.get("breadth", ""),
                    "mapped_reads": row.get("reads", ""),
                    "read_pairs_ge97": row.get("read_pairs_ge97", ""),
                    "assembly_support_bases": row.get("assembly_support_bases", ""),
                    "assembled_sequences": row.get("assembled_sequences", ""),
                    "assembled_sources": row.get("assembled_sources", ""),
                    "mag_from_library": row.get("mag_from_library", ""),
                }
            )
    context = []
    for sample, row in sorted(meta.items()):
        genus = row["species_current"].split()[0]
        tissue = row.get("tissue", "").strip()
        f = flow[sample]
        physical_count = int(f["flowcell_count"])
        eligible = f["inference_eligible"] == "true"
        if eligible != (physical_count == 1):
            raise ValueError(f"Single-flowcell eligibility conflict: {sample}")
        r = {
            "library_id": row["library_id"],
            "sample_id": sample,
            "species_current": row["species_current"],
            "host_genus": genus,
            "host_group": "Physalia" if genus == "Physalia" else "non_Physalia",
            "lifestyle": "neustonic"
            if genus == "Physalia"
            else "benthic"
            if row["species_current"] == "Stephalia dilata"
            else "pelagic",
            **{
                key: row.get(key, "")
                for key in (
                    "study",
                    "ocean_region",
                    "locality",
                    "collection_date",
                    "latitude",
                    "longitude",
                )
            },
            "collection_year": row.get("collection_date", "")[:4],
            "tissue_raw": tissue,
            "tissue_group": tissue.casefold() if tissue else "unknown",
            **{
                key: row.get(key, "")
                for key in (
                    "depth_m",
                    "depth_min_m",
                    "depth_max_m",
                    "depth_original",
                    "collection_depth_basis",
                    "collection_depth_source",
                )
            },
            **depth_context(row),
            "flowcell_count": physical_count,
            "flowcells": f["flowcells"],
            "single_flowcell_eligible": str(eligible).lower(),
            "read_pairs": row.get("read_pairs", ""),
        }
        for bucket, count_column in COUNTS.items():
            r[count_column] = len(present[sample][bucket])
            r[f"{bucket}_detected"] = int(bool(present[sample][bucket]))
        r["dt68_grade"] = (
            "high_confidence"
            if "high_confidence" in grades[sample]
            else "validated"
            if grades[sample]
            else "not_detected"
        )
        for genome in FOCAL[:-1]:
            r[f"{genome}_detected"] = int(bool(present[sample][genome]))
        context.append(r)
    return context, evidence


def summarize(rows: list[dict], group_columns: tuple[str, ...]) -> list[dict]:
    """Summarize positive libraries and target/group counts using all denominators."""
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[key] for key in group_columns)].append(row)
    output = []
    for keys, subset in sorted(groups.items()):
        result = dict(zip(group_columns, keys))
        result["total_libraries"] = len(subset)
        result["single_flowcell_libraries"] = sum(
            row["single_flowcell_eligible"] == "true" for row in subset
        )
        result["library_ids"] = ";".join(row["library_id"] for row in subset)
        for bucket, count_column in COUNTS.items():
            result[f"{bucket}_positive_libraries"] = sum(
                row[f"{bucket}_detected"] for row in subset
            )
            result[count_column] = sum(row[count_column] for row in subset)
        output.append(result)
    return output


def regional_incidence(rows: list[dict]) -> list[dict]:
    """Compare all Physalia and tentacle-only fractions without new inference."""
    physalia = [row for row in rows if row["host_group"] == "Physalia"]
    regions = ["all_regions"] + sorted(
        {row["ocean_region"] or "unknown" for row in physalia}
    )
    output = []
    for scope in ("all_libraries", "single_flowcell"):
        scoped = [
            row
            for row in physalia
            if scope == "all_libraries" or row["single_flowcell_eligible"] == "true"
        ]
        for tissue in ("all_tissues", "tentacle"):
            tissue_rows = [
                row
                for row in scoped
                if tissue == "all_tissues" or row["tissue_group"] == "tentacle"
            ]
            for region in regions:
                subset = [
                    row
                    for row in tissue_rows
                    if region == "all_regions"
                    or (row["ocean_region"] or "unknown") == region
                ]
                for target in FOCAL:
                    column = (
                        "vibrio_detected"
                        if target == "Vibrio"
                        else f"{target}_detected"
                    )
                    positives = [row for row in subset if row[column]]
                    output.append(
                        {
                            "scope": scope,
                            "tissue": tissue,
                            "ocean_region": region,
                            "target": target,
                            "positive_libraries": len(positives),
                            "total_libraries": len(subset),
                            "fraction_positive": f"{len(positives) / len(subset):.8g}"
                            if subset
                            else "",
                            "positive_library_ids": ";".join(
                                row["library_id"] for row in positives
                            ),
                            "denominator_library_ids": ";".join(
                                row["library_id"] for row in subset
                            ),
                        }
                    )
    return output


def sha256(path: Path) -> str:
    """Hash the small inputs and outputs supporting the retained comparison."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(root: Path, output: Path) -> dict:
    """Execute the bounded descriptive analysis and capture contemporaneous evidence."""
    start = datetime.now(timezone.utc)
    paths = {key: root / relative for key, relative in INPUTS.items()}
    inputs = {key: read_table(path) for key, path in paths.items()}
    context, evidence = derive_context(inputs)
    if len(context) != 205:
        raise ValueError(
            "Review descriptive scope if the accepted 205-library cohort changes"
        )
    output.mkdir(parents=True, exist_ok=True)
    tables = {
        "library_context.tsv": context,
        "associate_detections.tsv": evidence,
        "depth_summary.tsv": summarize(context, ("depth_group",)),
        "host_detection_summary.tsv": summarize(context, ("species_current",)),
        "tissue_summary.tsv": summarize(context, ("host_group", "tissue_group")),
        "dt68_specimens.tsv": [
            row for row in context if row["host_genus"] in {"Nanomia", "Resomia"}
        ],
        "dt68_depth_summary.tsv": summarize(
            [row for row in context if row["host_genus"] in {"Nanomia", "Resomia"}],
            ("species_current", "depth_group"),
        ),
        "tissue_region_incidence.tsv": regional_incidence(context),
        "tissue_region_host_incidence.tsv": [
            {"species_current": "Physalia utriculus", **row}
            for row in regional_incidence(
                [
                    row
                    for row in context
                    if row["species_current"] == "Physalia utriculus"
                ]
            )
        ],
    }
    for filename, rows in tables.items():
        write_table(output / filename, rows)
    summary = {
        "libraries": len(context),
        "single_flowcell_libraries": sum(
            row["single_flowcell_eligible"] == "true" for row in context
        ),
        "depth_status_counts": dict(
            sorted(Counter(row["depth_status"] for row in context).items())
        ),
        "depth_group_counts": dict(
            sorted(Counter(row["depth_group"] for row in context).items())
        ),
        "cohort": summarize(context, ())[0],
        "host_group_summary": summarize(context, ("host_group",)),
        "tissue_summary": tables["tissue_summary.tsv"],
        "dt68_depth_summary": tables["dt68_depth_summary.tsv"],
        "methods": {
            "scope": "All 205 libraries descriptive; additional single-physical-flowcell descriptive subsets only; no new tests or models.",
            "positive_grades": sorted(VALID),
            "eukaryote_selection": "count_as_detection=true in reviewed evidence; no ancestor/descendant double-counting.",
            "units": "Bacterial/archaeal/viral counts are accepted genome targets; eukaryote counts are reviewed reporting groups. They are not equivalent richness measures.",
            "zero_interpretation": "No accepted detection in the reviewed evidence, not established biological absence.",
            "depth_strata": "Physalia at assigned/recorded surface; other hosts entirely within 0-20 m, entirely below 20 m, crossing 20 m, or unknown. Intervals have no imputed midpoint.",
            "physalia_depths": "150 curator-assigned surface values and one source-record zero are distinguished.",
            "roles": "Roles are evidence annotations. A library may contain both potential prey and potential parasites.",
            "tissue_groups": "Source tissue descriptions are lowercased; empty values are unknown. No anatomical localization is inferred from DNA presence.",
            "host_tissue_check": "A fixed secondary descriptive subset of Physalia utriculus compares regional incidence within one host species and tissue description; all-library and single-flowcell scopes remain separate.",
            "lifestyle": "Physalia neustonic; Stephalia dilata benthic; remaining taxa broadly pelagic, not depth-derived.",
            "interpretation": "Post hoc exploratory descriptions; host, depth, region, tissue, sampling time and sequencing batch overlap.",
        },
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    source = Path(__file__).resolve()
    snapshot = output / "provenance/analyze_collection_context.py"
    snapshot.parent.mkdir(parents=True, exist_ok=True)
    snapshot.write_bytes(source.read_bytes())
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True
    ).strip()
    dirty = subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=root, text=True
    ).splitlines()

    def relative(path: Path) -> str:
        try:
            return str(path.resolve().relative_to(root.resolve()))
        except ValueError:
            return str(path.resolve())

    provenance = {
        "analysis_id": "collection-context-descriptive",
        "run_id": start.strftime("%Y%m%dT%H%M%S.%fZ"),
        "started_at_utc": start.isoformat(),
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "complete",
        "command": shlex.join([sys.executable, *sys.argv]),
        "working_directory": str(Path.cwd()),
        "python_version": platform.python_version(),
        "dependencies": "Python standard library only",
        "source_base_git_revision": revision,
        "source_dirty_state": dirty,
        "script": relative(source),
        "script_sha256": sha256(source),
        "source_snapshot": relative(snapshot),
        "inputs_sha256": {relative(path): sha256(path) for path in paths.values()},
        "outputs_sha256": {
            relative(output / name): sha256(output / name)
            for name in [*tables, "summary.json"]
        },
        "parameters": summary["methods"],
        "checks": {
            "one_to_one_identifiers": True,
            "all_libraries_retained": len(context),
            "all_evidence_samples_joined": True,
            "depth_intervals_valid": True,
            "flowcell_policy_verified": True,
            "evidence_rows": len(evidence),
        },
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return summary


def main() -> None:
    """CLI entry point; defaults are explicit accepted publication inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    output = args.output_dir or args.root / "figures/collection_context"
    result = run(args.root, output)
    print(
        json.dumps(
            {
                "libraries": result["libraries"],
                "depth_group_counts": result["depth_group_counts"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
