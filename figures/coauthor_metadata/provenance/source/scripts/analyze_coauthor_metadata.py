#!/usr/bin/env python3
"""Describe tissue/depth sampling and prepare bounded bacterial comparisons."""

from __future__ import annotations

import argparse
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from analyze_collection_context import (
    COUNTS,
    INPUTS,
    derive_context,
    read_table,
    summarize,
    write_table,
    sha256,
)


def join_mapping(rows: list[dict], mapping: list[dict]) -> list[dict]:
    """Require exactly one actual mapping-input denominator for each specimen."""
    by_sample = {row["sample_id"]: row for row in mapping}
    if len(by_sample) != len(mapping) or set(by_sample) != {
        row["sample_id"] for row in rows
    }:
        raise ValueError("Mapping summary and specimen metadata must match one-to-one")
    for row in rows:
        source = by_sample[row["sample_id"]]
        pairs = int(source["input_pairs"])
        if pairs <= 0:
            raise ValueError(f"No mapping input pairs: {row['sample_id']}")
        row.update(
            input_pairs=pairs,
            log10_input_pairs=math.log10(pairs),
            host_route=source["host_route"],
            read_source=source["read_source"],
            background_control_flag=row["sample_id"] == "Manko2026__T676-SS3",
        )
    return rows


def select(rows: list[dict], filters: dict[str, list[str]]) -> list[dict]:
    """Select declared biological subsets before considering detection outcomes."""
    return [
        row
        for row in rows
        if all(row[key] in values for key, values in filters.items())
    ]


def join_fish(rows: list[dict], fish: list[dict]) -> list[dict]:
    """Retain the separate marker screen's counts and sequencing denominator."""
    indexed = {row["analysis_sample_id"]: row for row in fish}
    if len(indexed) != len(fish) or set(indexed) != {row["sample_id"] for row in rows}:
        raise ValueError("Fish screen and specimen metadata must match one-to-one")
    for row in rows:
        source = indexed[row["sample_id"]]
        row["fish_12s_pairs_screened"] = int(source["pairs_screened"])
        row["fish_12s_matching_pairs"] = int(source["long_close_fish_match"])
        row["fish_12s_detected"] = int(row["fish_12s_matching_pairs"] > 0)
    return rows


def label_allocations(rows: list[dict], group: str) -> int:
    """Count distinct group-label allocations under physical-flowcell blocking."""
    blocks = defaultdict(Counter)
    for row in rows:
        blocks[row["flowcells"]][row[group]] += 1
    allocations = 1
    for counts in blocks.values():
        ways = math.factorial(sum(counts.values()))
        for count in counts.values():
            ways //= math.factorial(count)
        allocations *= ways
    return allocations


def group_description(rows: list[dict], columns: tuple[str, ...]) -> list[dict]:
    """Add continuous depth and actual sequencing effort to incidence summaries."""
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[column] for column in columns)].append(row)
    output = summarize(rows, columns)
    for result in output:
        members = grouped[tuple(result[column] for column in columns)]
        depths = sorted(float(row["depth_m"]) for row in members if row["depth_m"])
        result.update(
            numeric_depths=len(depths),
            depth_min=min(depths) if depths else "",
            depth_max=max(depths) if depths else "",
            input_pairs_min=min(row["input_pairs"] for row in members),
            input_pairs_max=max(row["input_pairs"] for row in members),
            fish_12s_positive_libraries=sum(
                row["fish_12s_detected"] for row in members
            ),
            fish_12s_matching_pairs=sum(
                row["fish_12s_matching_pairs"] for row in members
            ),
            physical_flowcells=len(
                {
                    row["flowcells"]
                    for row in members
                    if row["single_flowcell_eligible"] == "true"
                }
            ),
        )
    return output


def depth_design(rows: list[dict]) -> list[dict]:
    """Expose within-host replication and depth/locality overlap without a model."""
    output = []
    by_host = defaultdict(list)
    for row in rows:
        by_host[row["species_current"]].append(row)
    for host, members in sorted(by_host.items()):
        eligible = [row for row in members if row["single_flowcell_eligible"] == "true"]
        depths = {row["depth_m"] for row in members if row["depth_m"]}
        groups_by_locality = defaultdict(set)
        for row in members:
            groups_by_locality[row["locality"]].add(row["depth_group"])
        crosses = sum(len(values) > 1 for values in groups_by_locality.values())
        if host.startswith("Physalia "):
            status = "no_collection_depth_variation"
        elif len(members) < 4 or len(eligible) < 4:
            status = "insufficient_within_host_single_flowcell_replication"
        elif len(depths) < 2:
            status = "no_collection_depth_variation"
        else:
            status = "heterogeneous_host_identification_or_unresolved_confounding"
        if host == "Nanomia septata":
            status = "shallow_deep_contrast_confounded_with_locality"
        output.append(
            {
                "species_current": host,
                "libraries": len(members),
                "single_flowcell_libraries": len(eligible),
                "distinct_point_or_derived_depths": len(depths),
                "localities_spanning_depth_groups": crosses,
                "depth_min_m": min(float(d) for d in depths) if depths else "",
                "depth_max_m": max(float(d) for d in depths) if depths else "",
                "status": status,
            }
        )
    return output


def summarize_host_references(
    config: dict, context: list[dict]
) -> tuple[list[dict], list[dict]]:
    """Keep current read-removal summaries distinct from historical v2 controls."""
    metadata = {row["sample_id"]: row for row in context}
    stats = read_table(Path(config["host_depletion_stats"]))
    provenance = json.loads(Path(config["host_depletion_provenance"]).read_text())
    if provenance["status"] != "complete" or provenance["outputs_sha256"][
        config["host_depletion_stats"]
    ] != sha256(Path(config["host_depletion_stats"])):
        raise ValueError(
            "Host-depletion statistics do not match completed source provenance"
        )
    expected = {
        sample for sample, row in metadata.items() if row["host_route"] != "none"
    }
    if len(stats) != len(expected) or {row["sample"] for row in stats} != expected:
        raise ValueError(
            "Host-depletion summary must cover all and only reference-mapped libraries"
        )
    grouped = defaultdict(list)
    for row in stats:
        source = metadata[row["sample"]]
        ingested, remaining = int(row["reads_ingested"]), int(row["reads_remaining"])
        removed = 100 * (ingested - remaining) / ingested
        if (
            not 0 <= remaining <= ingested
            or abs(removed - float(row["percent_removed"])) > 0.000051
        ):
            raise ValueError(f"Inconsistent host-depletion counts: {row['sample']}")
        values = (removed, ingested, remaining)
        grouped[("species", source["species_current"])].append(values)
        group = (
            "Physalia"
            if source["host_group"] == "Physalia"
            else source["species_current"]
        )
        grouped[("reference_group", group)].append(values)
    summary = []
    for (scope, host), members in sorted(grouped.items()):
        percentages = [row[0] for row in members]
        summary.append(
            {
                "scope": scope,
                "host": host,
                "libraries": len(members),
                "percent_removed_min": min(percentages),
                "percent_removed_median": statistics.median(percentages),
                "percent_removed_max": max(percentages),
                "reads_ingested": sum(row[1] for row in members),
                "reads_remaining": sum(row[2] for row in members),
            }
        )
    historical = read_table(Path(config["historical_host_concordance"]))
    keys = [(row["sample_id"], row["domain"]) for row in historical]
    if len(keys) != len(set(keys)) or len(keys) != 20:
        raise ValueError(
            "Expected two domains in the historical ten-library comparison"
        )
    historical_groups = defaultdict(list)
    for row in historical:
        sample, domain = row["sample_id"], row["domain"]
        both, depleted, full = (
            int(row[k]) for k in ("both", "only_host_depleted", "only_full_reads")
        )
        if (
            int(row["present_host_depleted"]) != both + depleted
            or int(row["present_full_reads"]) != both + full
        ):
            raise ValueError(f"Historical concordance does not reconcile: {sample}")
        group = (
            "Physalia"
            if metadata[sample]["host_group"] == "Physalia"
            else metadata[sample]["species_current"]
        )
        historical_groups[("all", domain)].append(row)
        historical_groups[(group, domain)].append(row)
    paired = []
    for (host, domain), members in sorted(historical_groups.items()):
        paired.append(
            {
                "host_group": host,
                "domain": domain,
                "libraries": len(members),
                "catalog_context": "historical_fixed_v2",
                "target_scope": "combined_associate_and_endogenous_candidate_viral_catalog"
                if domain == "viral"
                else "bacterial_catalog_after_reviewed_domain_exclusions",
                **{
                    key: sum(int(row[key]) for row in members)
                    for key in (
                        "present_host_depleted",
                        "present_full_reads",
                        "both",
                        "only_host_depleted",
                        "only_full_reads",
                    )
                },
                "sample_ids": ";".join(sorted(row["sample_id"] for row in members)),
            }
        )
    return summary, paired


def run(config_path: Path) -> None:
    """Read completed mapping/evidence tables; never open sequence inputs."""
    config = json.loads(config_path.read_text())
    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    inputs = {key: read_table(Path(path)) for key, path in INPUTS.items()}
    context, evidence = derive_context(inputs)
    context = join_mapping(context, read_table(Path(config["mapping_summary"])))
    context = join_fish(context, read_table(Path(config["fish_summary"])))
    if len(context) != config["expected_libraries"]:
        raise ValueError("Unexpected cohort size")
    if any(not row["depth_m"] for row in context):
        raise ValueError("Expected the reviewed, depth-complete metadata table")
    tables = {
        "specimen_context.tsv": context,
        "associate_detections.tsv": evidence,
        "tissue_incidence.tsv": group_description(
            context, ("host_group", "tissue_group")
        ),
        "depth_incidence.tsv": group_description(context, ("study", "depth_group")),
        "host_depth_design.tsv": depth_design(context),
        "host_locality_incidence.tsv": group_description(
            context, ("species_current", "locality")
        ),
    }
    regional = []
    for scope in ("all_libraries", "single_flowcell"):
        scoped = [
            row
            for row in context
            if scope == "all_libraries" or row["single_flowcell_eligible"] == "true"
        ]
        for tissue in ("all_tissues", "tentacle"):
            selected = [
                row
                for row in scoped
                if tissue == "all_tissues" or row["tissue_group"] == tissue
            ]
            for subset_id, filters, groups in [
                (
                    "utriculus_regions",
                    {"species_current": ["Physalia utriculus"]},
                    ("ocean_region",),
                ),
                (
                    "sw_pacific_hosts",
                    {"host_group": ["Physalia"], "ocean_region": ["SW Pacific"]},
                    ("species_current",),
                ),
            ]:
                for result in group_description(select(selected, filters), groups):
                    regional.append(
                        {
                            "subset": subset_id,
                            "scope": scope,
                            "tissue": tissue,
                            "ocean_region": "",
                            "species_current": "",
                            **result,
                        }
                    )
    tables["targeted_incidence.tsv"] = regional
    designs, memberships = [], []
    for spec in config["contrasts"]:
        all_rows = select(context, spec["filters"])
        eligible = [
            row for row in all_rows if row["single_flowcell_eligible"] == "true"
        ]
        allocations = label_allocations(eligible, spec["group"])
        designs.append(
            {
                "contrast": spec["id"],
                "group_column": spec["group"],
                "all_libraries": len(all_rows),
                "single_flowcell_libraries": len(eligible),
                "distinct_label_allocations": str(allocations),
                "nominal_5pct_test_possible": allocations >= 20,
            }
        )
        for row in all_rows:
            memberships.append(
                {
                    "contrast": spec["id"],
                    "group": row[spec["group"]],
                    "sample_id": row["sample_id"],
                    "flowcell": row["flowcells"],
                    "single_flowcell_eligible": row["single_flowcell_eligible"],
                    "input_pairs": row["input_pairs"],
                }
            )
    tables["contrast_design.tsv"] = designs
    tables["contrast_membership.tsv"] = memberships
    host_summary, historical_pairs = summarize_host_references(config, context)
    tables["host_depletion_summary.tsv"] = host_summary
    tables["historical_host_concordance.tsv"] = historical_pairs
    (output / "host_reference_values.json").write_text(
        json.dumps(
            {
                "current_host_depletion": host_summary,
                "historical_fixed_v2_concordance": historical_pairs,
                "interpretation": "Removal percentages describe current read processing; paired control counts reuse the historical fixed v2 catalog and are not v3 results. Viral paired controls include both associate and endogenous-candidate targets. Neither comparison isolates a causal effect on microbial discovery.",
            },
            indent=2,
        )
        + "\n"
    )
    for name, rows in tables.items():
        write_table(output / name, rows)
    values = {
        "libraries": len(context),
        "numeric_depths": sum(bool(row["depth_m"]) for row in context),
        "tissue_unknown": sum(row["tissue_group"] == "unknown" for row in context),
        "fish_12s_positive_libraries": sum(row["fish_12s_detected"] for row in context),
        "fish_12s_matching_pairs": sum(
            row["fish_12s_matching_pairs"] for row in context
        ),
        "tissue_counts": {
            str(key): count
            for key, count in sorted(
                Counter(
                    (row["host_group"], row["tissue_group"]) for row in context
                ).items()
            )
        },
        "depth_design": tables["host_depth_design.tsv"],
        "contrast_design": designs,
        "detection_columns": COUNTS,
        "interpretation": config["interpretation"],
    }
    (output / "descriptive_values.json").write_text(json.dumps(values, indent=2) + "\n")
    print(json.dumps({"libraries": len(context), "contrast_design": designs}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, default=Path("config/coauthor_metadata_analysis.json")
    )
    run(parser.parse_args().config)
