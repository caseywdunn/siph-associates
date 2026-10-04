#!/usr/bin/env python3
"""Test physical-flowcell concentration within host-species x region strata.

Only libraries eligible in the normalized metadata enter inference. Descriptive
presence counts retain all libraries. Strata without two physical flowcells do
not enter the concentration statistic; a constant permutation distribution is
reported as untestable. Permutations and BH adjustment are deterministic.
"""

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

PRESENT = {"validated", "high_confidence"}
FIELDS = [
    "target_id", "taxonomy", "validated_presences", "eligible_presences",
    "permutable_presences", "excluded_presences", "flowcells",
    "observed_concentration", "p_value", "q_value", "status",
    "probable_contaminant", "identity_annotation", "analysis_subset",
    "grade_column",
]


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def concentration(labels: dict[str, str], members: set[str]) -> float:
    counts = Counter(labels[sample] for sample in members)
    return sum(count * count for count in counts.values()) / len(members)


def contamination_results(
    rule: dict, grades: list[dict[str, str]], metadata: list[dict[str, str]],
    grade_column: str = "grade", subset: str = "all",
) -> list[dict]:
    """Return separate descriptive and inferential counts for each tested genome."""
    meta = {
        row["sample_id"]: row for row in metadata
        if subset == "all" or row["host_species"].startswith("Nanomia")
    }
    libraries = {
        sample: row for sample, row in meta.items()
        if row["inference_eligible"] == "true"
    }
    strata = defaultdict(list)
    for sample, row in sorted(libraries.items()):
        strata[(row["host_species"], row["ocean_region"] or "unknown")].append(sample)
    permutable = {
        key: members for key, members in sorted(strata.items())
        if len({libraries[sample]["flowcell"] for sample in members})
        >= rule["min_flowcells_in_stratum"]
    }
    permutable_samples = {sample for members in permutable.values() for sample in members}
    present, taxonomy = defaultdict(set), {}
    for row in grades:
        if (row["sample_id"] in meta and row["role"] != "decoy"
                and row[grade_column] in PRESENT):
            present[row["target_id"]].add(row["sample_id"])
            taxonomy[row["target_id"]] = row["taxonomy"]

    generator = random.Random(rule["seed"])
    base = {sample: libraries[sample]["flowcell"] for sample in sorted(permutable_samples)}
    permutations = []
    for _ in range(rule["permutations"]):
        labels = dict(base)
        for members in permutable.values():
            shuffled = [base[sample] for sample in members]
            generator.shuffle(shuffled)
            labels.update(zip(members, shuffled))
        permutations.append(labels)

    results = []
    for target, all_members in sorted(present.items()):
        if len(all_members) < rule["min_validated_presences"]:
            continue
        eligible_members = all_members & libraries.keys()
        members = all_members & permutable_samples
        observed = concentration(base, members) if members else None
        status, p_value = "untestable_no_permutable_presences", None
        if len(eligible_members) < rule["min_validated_presences"]:
            status = "below_minimum_eligible_presences"
        elif members:
            null = [concentration(labels, members) for labels in permutations]
            # The statistic may be fixed even in a stratum with multiple flowcells.
            # For example, all its libraries can carry the genome.
            if min(null + [observed]) == max(null + [observed]):
                status = "untestable_constant_statistic"
            else:
                status = "tested"
                p_value = (sum(value >= observed for value in null) + 1) / (len(null) + 1)
        genus = next((rank[3:] for rank in taxonomy[target].split(";")
                      if rank.startswith("g__")), "")
        results.append({
            "target_id": target, "taxonomy": taxonomy[target],
            "validated_presences": len(all_members),
            "eligible_presences": len(eligible_members),
            "permutable_presences": len(members),
            "excluded_presences": len(all_members - eligible_members),
            "flowcells": len({libraries[sample]["flowcell"] for sample in members}),
            "observed_concentration": "" if observed is None else observed,
            "p_value": "" if p_value is None else p_value,
            "status": status,
            "identity_annotation": str(genus.split("_")[0]
                                       in rule["identity_annotation_genera"]).lower(),
            "analysis_subset": subset, "grade_column": grade_column,
        })
    tested = sorted((row for row in results if row["status"] == "tested"),
                    key=lambda row: row["p_value"])
    running = 1.0
    for rank, row in reversed(list(enumerate(tested, start=1))):
        running = min(running, row["p_value"] * len(tested) / rank)
        row["q_value"] = running
    for row in results:
        row.setdefault("q_value", "")
        row["probable_contaminant"] = str(
            row["status"] == "tested" and row["q_value"] < rule["fdr"]
        ).lower()
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", required=True, type=Path)
    parser.add_argument("--grades", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--grade-column", default="grade",
                        choices=("grade", "grade_at_5pct", "grade_at_20pct"))
    parser.add_argument("--subset", default="all", choices=("all", "nanomia"))
    args = parser.parse_args()
    rule = json.loads(args.analysis.read_text())["contamination_test"]
    results = contamination_results(rule, rows(args.grades), rows(args.metadata),
                                 args.grade_column, args.subset)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    temporary.replace(args.output)
    print(f"genomes considered={len(results)} "
          f"tested={sum(row['status'] == 'tested' for row in results)} "
          f"flagged={sum(row['probable_contaminant'] == 'true' for row in results)}")


if __name__ == "__main__":
    main()
