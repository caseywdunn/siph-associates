"""Stage aggregate human-screen results without candidate sequences or alignments."""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import statistics
import shlex
import shutil
import sys
from pathlib import Path


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stage(config_path, confirmation_scope, execution_ids, output):
    config = json.loads(config_path.read_text())
    root = Path(config["run"])
    sample_path = root / "cohort_summary.tsv"
    control_path = root / "controls.tsv"
    controls = json.loads((root / "controls.passed.json").read_text())
    samples = read_table(sample_path)
    manifest_path = Path(config["manifest"])
    manifest = {row["sample_id"]: row for row in read_table(manifest_path)}
    if len(samples) != len(manifest) or {row["sample_id"] for row in samples} != set(
        manifest
    ):
        raise ValueError("Cohort report does not cover the source inventory exactly")
    if controls["human_recovery"] < 0.7 or controls["nonhuman_false_positives"]:
        raise ValueError("Controls have not passed")
    confirmation = root / "confirmation" / confirmation_scope
    nt = json.loads((confirmation / "taxonomic_review.json").read_text())
    selection = json.loads((confirmation / "selection.json").read_text())
    executions = []
    for job_id in execution_ids:
        directory = root.parent / "executions" / job_id
        if (directory / "exit_code.txt").read_text().strip() != "0":
            raise ValueError(f"Execution {job_id} has not completed successfully")
        execution = {
            "slurm_job_id": job_id,
            "source_sha256": sha256(directory / "source.tar.gz"),
            "started": (directory / "started.txt").read_text().strip(),
            "finished": (directory / "finished.txt").read_text().strip(),
        }
        if (directory / "reference_identity.json").exists():
            execution["reference_identity"] = json.loads(
                (directory / "reference_identity.json").read_text()
            )
        if (directory / "nt_database_info.txt").exists():
            execution["nt_database_info"] = (
                directory / "nt_database_info.txt"
            ).read_text()
        executions.append(execution)
    output.mkdir(parents=True, exist_ok=True)
    public_samples = []
    for row in samples:
        sample = row["sample_id"]
        source = manifest[sample]
        pair_counts = nt["per_sample_pair_classifications"].get(sample, {})
        assessed_pairs = sum(pair_counts.values())
        status = (
            "selected_candidates_reviewed_against_nt; mapping_count_not_taxonomically_adjusted"
            if assessed_pairs
            else "no_retained_candidates"
            if int(row["stringent_human_compatible_pairs"]) == 0
            else "not_in_nt_confirmation_scope"
        )
        public_samples.append(
            {
                "sample_id": sample,
                "study": source["study"],
                **{key: value for key, value in row.items() if key != "sample_id"},
                "taxonomic_confirmation": status,
                "distinct_candidate_pairs_reviewed_nt": assessed_pairs,
                "nt_both_mates_human_only_near_best": pair_counts.get(
                    "both_mates_human_only_near_best", 0
                ),
                "nt_human_with_alternatives": pair_counts.get(
                    "human_with_alternatives", 0
                ),
                "nt_incomplete_confirmation": pair_counts.get(
                    "incomplete_confirmation", 0
                ),
                "nt_unresolved_taxonomy_or_hit_limit": pair_counts.get(
                    "unresolved_truncation_or_unknown_taxonomy", 0
                ),
                "nt_conflicting_mates": pair_counts.get("conflicting_mates", 0),
                "nt_human_not_supported": pair_counts.get(
                    "human_not_supported_by_nt", 0
                ),
            }
        )
    with (output / "sample_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(public_samples[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(public_samples)
    control_rows = read_table(control_path)
    for row in control_rows:
        row["taxonomic_confirmation"] = "known_origin_simulated_control"
    with (output / "control_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(control_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(control_rows)
    # Taxonomic assets contain only aggregate categories, never human sequences,
    # candidate identifiers, genomic coordinates or individual-level information.
    (output / "nt_summary.json").write_text(
        json.dumps({"scope": confirmation_scope, **selection, **nt}, indent=2) + "\n"
    )
    count_keys = [
        "raw_selected_pairs",
        "qc_pairs",
        "bowtie2_concordant_pairs",
        "stringent_human_compatible_pairs",
        "distinct_sequence_pairs",
    ]
    totals = {key: sum(int(row[key]) for row in samples) for key in count_keys}
    maximum = max(samples, key=lambda r: float(r["percent_selected_raw_pairs"] or 0))
    formalin = next(row for row in samples if row["sample_id"] == "Manko2026__T676-SS3")
    summary = {
        "sampling_design": "equal per-source-file prefixes; not an unbiased whole-library abundance estimate",
        "maximum_selected_raw_pairs_per_library": config["pairs_per_library"],
        "libraries": len(samples),
        "studies": dict(
            Counter(manifest[row["sample_id"]]["study"] for row in samples)
        ),
        **totals,
        "libraries_with_human_compatible_pairs": sum(
            int(row["stringent_human_compatible_pairs"]) > 0 for row in samples
        ),
        "libraries_without_qc_pairs": sum(int(row["qc_pairs"]) == 0 for row in samples),
        "highest_selected_prefix_fraction": {
            key: maximum[key]
            for key in (
                "sample_id",
                "stringent_human_compatible_pairs",
                "percent_selected_raw_pairs",
                "occupied_1Mb_windows",
            )
        },
        "T676_SS3": {
            key: formalin[key]
            for key in count_keys
            + [
                "percent_selected_raw_pairs",
                "occupied_1Mb_windows",
                "reference_sequences",
            ]
        },
        "controls": controls,
        "taxonomic_confirmation_scope": confirmation_scope,
        "taxonomic_confirmation": nt,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    nt_counts = nt["pair_classification_counts"]
    macros = {
        "currentHumanNtPairs": f"{nt['selected_pair_occurrences']:,}",
        "currentHumanNtHumanOnlyPairs": f"{nt_counts.get('both_mates_human_only_near_best', 0):,}",
        "currentHumanNtAlternativePairs": f"{nt_counts.get('human_with_alternatives', 0):,}",
        "currentHumanNtIncompletePairs": f"{nt_counts.get('incomplete_confirmation', 0):,}",
        "currentHumanNtCappedOrUnknownPairs": f"{nt_counts.get('unresolved_truncation_or_unknown_taxonomy', 0):,}",
        "currentHumanNtConflictingPairs": f"{nt_counts.get('conflicting_mates', 0):,}",
        "currentHumanNtNonhumanPairs": f"{nt_counts.get('human_not_supported_by_nt', 0):,}",
        "currentHumanMaximumPairs": f"{int(maximum['stringent_human_compatible_pairs']):,}",
        "currentHumanMaximumPercent": f"{float(maximum['percent_selected_raw_pairs']):.4f}",
        "currentHumanNegativeControlPairs": f"{controls['negative_pairs']:,}",
        "currentHumanScreenLibraries": len(samples),
        "currentHumanScreenPairs": f"{totals['raw_selected_pairs']:,}",
        "currentHumanConcordantPairs": f"{totals['bowtie2_concordant_pairs']:,}",
        "currentHumanMedianConcordantPairs": f"{statistics.median(int(row['bowtie2_concordant_pairs']) for row in samples):g}",
        "currentHumanMedianStringentPairs": f"{statistics.median(int(row['stringent_human_compatible_pairs']) for row in samples):g}",
        "currentHumanMedianPairs": f"{statistics.median(int(row['stringent_human_compatible_pairs']) for row in samples):g}",
        "currentHumanCompatiblePairs": f"{totals['stringent_human_compatible_pairs']:,}",
        "currentHumanPositiveLibraries": summary[
            "libraries_with_human_compatible_pairs"
        ],
        "currentHumanFormalinPairs": formalin["stringent_human_compatible_pairs"],
        "currentHumanFormalinPercent": f"{float(formalin['percent_selected_raw_pairs']):.4f}",
        "currentHumanControlRecoveryPercent": f"{100 * controls['human_recovery']:.2f}",
    }
    (output / "values.tex").write_text(
        "% Human-compatible mapping metrics in source-file prefixes.\n"
        + "".join(
            f"\\newcommand{{\\{key}}}{{{value}}}\n" for key, value in macros.items()
        )
    )
    inputs = [
        Path(__file__),
        Path("scripts/human_screen/review_nt.py"),
        config_path,
        manifest_path,
        sample_path,
        control_path,
        root / "controls.passed.json",
        confirmation / "selection.json",
        confirmation / "taxonomic_review.json",
    ]
    source_files = [
        config_path,
        Path("envs/human-screen.yaml"),
        Path("workflow/human_screen.smk"),
        Path("workflow/human_confirmation.smk"),
        *sorted(Path("scripts/human_screen").glob("*.py")),
        *sorted(Path("scripts/human_screen").glob("*.sbatch")),
    ]
    for path in source_files:
        destination = output / "provenance/source" / path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, destination)
    provenance = {
        "analysis_id": "human-compatible-prefix-audit",
        "run_id": datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "status": "complete",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": shlex.join([sys.executable, *sys.argv]),
        "report_runtime": {
            "python": sys.version,
            "snakemake": __import__("snakemake").__version__,
            "pysam": __import__("pysam").__version__,
            "biopython": __import__("Bio").__version__,
        },
        "source_base_git_revision": (
            root.parent / "executions" / execution_ids[0] / "git_head.txt"
        )
        .read_text()
        .strip(),
        "source_snapshot": str(output / "provenance/source"),
        "source_sha256": {str(path): sha256(path) for path in source_files},
        "inputs_sha256": {
            str(path): sha256(path) for path in inputs if path not in source_files
        },
        "outputs_sha256": {
            str(path): sha256(path)
            for path in sorted(output.iterdir())
            if path.is_file() and path.name != "provenance.json"
        },
        "executions": executions,
        "inputs": [{"path": str(path), "sha256": sha256(path)} for path in inputs],
        "outputs": [
            {"path": str(path.name), "sha256": sha256(path)}
            for path in sorted(output.iterdir())
            if path.is_file() and path.name != "provenance.json"
        ],
    }
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument(
        "--confirmation-scope", choices=("pilot", "cohort"), required=True
    )
    parser.add_argument("--execution-ids", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    stage(args.config, args.confirmation_scope, args.execution_ids, args.output)


if __name__ == "__main__":
    main()
