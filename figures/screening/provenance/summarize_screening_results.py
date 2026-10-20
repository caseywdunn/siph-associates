#!/usr/bin/env python3
"""Summarize retained screening and validation outputs for the manuscript.

This reporting command reads existing tables only. It does not run classifiers,
align reads, change presence grades, or refresh other publication exports.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/screening"
INPUTS: dict[str, str] = {}
ACCEPTED = {"validated", "high_confidence"}


def sha256(path: Path) -> str:
    """Hash files incrementally, including the large nomination table."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source(relative: str) -> Path:
    path = ROOT / relative
    INPUTS[relative] = sha256(path)
    return path


def read(relative: str) -> list[dict[str, str]]:
    with source(relative).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name: str, rows: list[dict]) -> None:
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def controls(relative: str, identity: int) -> dict[str, Counter]:
    """Count the same control lineages before and after competitive remapping."""
    result = {name: Counter() for name in ("Insecta", "Embryophyta", "Mammalia")}
    for row in read(relative):
        if int(row["min_identity_band"]) < identity:
            continue
        ranks = row["lineage"].split(";")
        if "Homo" in ranks:
            continue
        for name, counts in result.items():
            if name in ranks:
                counts[row["sample_id"]] += int(row["read_pairs"])
    return result


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    qc = read("data/results/aggregation/cohort/library_qc.tsv")
    samples = {row["sample_id"] for row in qc}
    if len(samples) != len(qc) or any(row["status"] != "PASS" for row in qc):
        raise ValueError("Screen scope must contain unique, accepted libraries")
    gate = json.loads(source("config/eukaryote_gate.json").read_text())
    identity = gate["presence"]["validated_reads"]["min_identity_band"]
    read_threshold = gate["presence"]["validated_reads"]["min_read_pairs"]
    analysis = json.loads(source("config/phase6_analysis.json").read_text())
    assembly_bases = analysis["evidence_grades"]["high_confidence_assembly_support"][
        "min_aligned_bases"
    ]

    counts: Counter = Counter()
    identifiers: dict[tuple[str, str], set[str]] = defaultdict(set)
    libraries: dict[tuple[str, str], set[str]] = defaultdict(set)
    sylph = []
    myco_labels = {method: Counter() for method in ("Bracken", "phyloFlash")}
    with source("data/results/aggregation/cohort/candidate_nominations.tsv").open(
        newline=""
    ) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sample_id"] not in samples:
                raise ValueError("Nomination outside accepted library scope")
            key = (row["source"], row["rank"])
            counts[key] += 1
            identifiers[key].add(row["candidate_id"])
            libraries[key].add(row["sample_id"])
            if row["source"] == "sylph":
                sylph.append(row)
            # Broad labels are evidence of a group signal, not identification of
            # the particular MAGs or a new organism-presence criterion.
            label = (row["candidate_id"] + " " + row["candidate_name"]).lower()
            if "mycoplasm" in label:
                if key == ("Bracken", "genus") and float(row["estimated_reads"]) >= 10:
                    myco_labels["Bracken"][row["sample_id"]] += 1
                elif row["source"] == "phyloFlash" and float(row["raw_support"]) >= 1:
                    myco_labels["phyloFlash"][row["sample_id"]] += 1
    method_rows = [
        {
            "method": method,
            "rank": rank,
            "assignment_rows": counts[(method, rank)],
            "distinct_labels": len(identifiers[(method, rank)]),
            "libraries": len(libraries[(method, rank)]),
        }
        for method, rank in sorted(counts)
    ]
    for rank, qc_field in (
        ("genus", "bracken_genus_rows"),
        ("species", "bracken_species_rows"),
    ):
        if counts[("Bracken", rank)] != sum(int(row[qc_field]) for row in qc):
            raise ValueError("Bracken nomination counts differ from accepted QC")
    nominations = read(
        "data/results/phase4_catalog/nomination_evidence/sylph_nominations.tsv"
    )
    for group in nominations:
        hits = [
            row
            for row in sylph
            if group["source"] == "all"
            or ("OceanDNA" if row["candidate_name"].startswith("OceanDNA") else "GTDB")
            == group["source"]
        ]
        observed = (
            len(hits),
            len({r["candidate_id"] for r in hits}),
            len({r["sample_id"] for r in hits}),
        )
        if observed != tuple(int(group[k]) for k in ("hits", "genomes", "libraries")):
            raise ValueError("Retained sylph summary disagrees with raw nominations")
    corroboration = read(
        "data/results/phase4_catalog/nomination_evidence/gtdb_corroboration.tsv"
    )
    uncorroborated = read(
        "data/results/phase4_catalog/nomination_evidence/gtdb_uncorroborated.tsv"
    )

    catalog = read("figures/catalog_review/bacterial_catalog.manifest.tsv")
    representatives = {
        row["catalog_id"]: row
        for row in catalog
        if row["status"] == "representative" and row["role"] != "decoy"
    }
    grades = [
        row
        for row in read("data/results/phase6_analysis/v2/grades/bacterial_grades.tsv")
        if row["role"] != "decoy"
    ]
    keys = {(row["sample_id"], row["target_id"]) for row in grades}
    if keys != {(s, t) for s in samples for t in representatives} or len(keys) != len(
        grades
    ):
        raise ValueError("Genome grades must cover each library/catalog pair once")
    nominated = [row for row in grades if row["nominated"] == "true"]
    nominated_targets = {row["target_id"] for row in nominated}
    mag_only = set(representatives) - nominated_targets
    if any(representatives[target]["role"] != "mag" for target in mag_only):
        raise ValueError("Non-nominated catalog targets are not all MAG-derived")
    accepted = [row for row in grades if row["grade"] in ACCEPTED]
    overlap = []
    for row in accepted:
        target = row["target_id"]
        group = (
            "same_library_sylph"
            if row["nominated"] == "true"
            else "other_library_sylph"
            if target in nominated_targets
            else "MAG_only"
        )
        overlap.append(
            {
                **row,
                "screen_overlap": group,
                "catalog_genome_id": representatives[target]["genome_id"],
                "cluster_members": representatives[target]["cluster_members"],
            }
        )
    myco = [row for row in overlap if "f__Metamycoplasmataceae;" in row["taxonomy"]]
    archaea = [row for row in overlap if row["taxonomy"].startswith("d__Archaea;")]
    if any(row["screen_overlap"] != "MAG_only" for row in myco + archaea):
        raise ValueError("Focal mycoplasmatacean/archaeal genomes are not all MAG-only")
    myco_samples = {row["sample_id"] for row in myco}
    if any(not myco_samples <= set(labels) for labels in myco_labels.values()):
        raise ValueError(
            "Broad Mycoplasma-related support is missing in a focal library"
        )
    myco_support = [
        {
            "sample_id": sample,
            "accepted_MAG_detections": sum(row["sample_id"] == sample for row in myco),
            "Bracken_related_genus_labels_ge10_reads": myco_labels["Bracken"][sample],
            "phyloFlash_related_lineage_labels_ge1_read": myco_labels["phyloFlash"][
                sample
            ],
        }
        for sample in sorted(myco_samples)
    ]
    mag_presences = sum(row["target_id"] in mag_only for row in accepted)
    dependence = read(
        "data/results/phase6_analysis/v2/sensitivity/leave_one_study_out.tsv"
    )
    prior = next(
        row
        for row in dependence
        if row["study"] == "MAG_only_genomes (no sylph nomination)"
    )
    if (len(mag_only), mag_presences) != (
        int(prior["genomes_nominated_only_by_study"]),
        int(prior["dependent_presences"]),
    ):
        raise ValueError(
            "MAG-only counts disagree with retained nomination-dependence summary"
        )

    before = controls("data/results/eukaryote_gate/read_evidence.tsv", identity)
    after = controls("data/results/eukaryote_gate/read_lineages.tsv", identity)
    calibration = read("data/results/eukaryote_gate/calibration/negative_controls.tsv")
    control_rows = []
    for name in before:
        for stage, group in (
            ("initial_phyloflash", before),
            ("competitive_host_aware", after),
        ):
            counts_by_library = group[name]
            if not set(counts_by_library) <= samples:
                raise ValueError("Control evidence outside accepted library scope")
            control_rows.append(
                {
                    "stage": stage,
                    "control": name,
                    "min_identity_percent": identity,
                    "libraries_with_pairs": sum(
                        n > 0 for n in counts_by_library.values()
                    ),
                    "max_pairs": max(counts_by_library.values(), default=0),
                    "libraries_at_presence_threshold": sum(
                        n >= read_threshold for n in counts_by_library.values()
                    ),
                }
            )
        retained = next(
            row
            for row in calibration
            if row["identity"] == f">={identity}%" and row["control"] == name
        )
        if (
            int(retained["libraries_with_pairs"]),
            int(retained["max_pairs_in_a_library"]),
        ) != (control_rows[-1]["libraries_with_pairs"], control_rows[-1]["max_pairs"]):
            raise ValueError("Competitive controls disagree with retained calibration")
    euk = read("data/results/eukaryote_gate/eukaryote_grades.tsv")
    euk_accepted = [row for row in euk if row["count_as_detection"] == "true"]
    euk_phyloflash = [
        row for row in euk_accepted if "phyloflash" in row["assembled_sources"]
    ]
    euk_shared = [
        row for row in euk_phyloflash if "assembly18s" in row["assembled_sources"]
    ]
    euk_assembly_only = [
        row for row in euk_accepted if row["assembled_sources"] == "assembly18s"
    ]
    rescue = [
        row
        for row in euk_assembly_only
        if row["sample_id"] == "Church2025__YPM-IZ-111009"
        and "Kinetoplastea" in row["reporting_unit"]
    ]
    if len(rescue) != 1 or int(rescue[0]["read_pairs_ge97"]) >= read_threshold:
        raise ValueError("Kinetoplastid assembly-rescue example no longer holds")
    totals = {row["source"]: row for row in nominations}
    corroborated = {
        row["evidence"]: int(row["hits_corroborated"]) for row in corroboration
    }
    values = {
        "Libraries": len(samples),
        "BrackenGenera": len(identifiers[("Bracken", "genus")]),
        "BrackenSpecies": len(identifiers[("Bracken", "species")]),
        "SylphGenomes": int(totals["all"]["genomes"]),
        "SylphLibraries": int(totals["all"]["libraries"]),
        "SylphHits": int(totals["all"]["hits"]),
        "GtdbGenomes": int(totals["GTDB"]["genomes"]),
        "OceanGenomes": int(totals["OceanDNA"]["genomes"]),
        "GtdbHits": int(totals["GTDB"]["hits"]),
        "BrackenCorroborated": corroborated["Bracken genus >=10 reads"],
        "PhyloflashCorroborated": corroborated["phyloFlash genus >=1 read"],
        "EitherCorroborated": corroborated["either"],
        "UncorroboratedGenomes": int(uncorroborated[0]["genomes"]),
        "NominatedPairs": len(nominated),
        "NominatedAccepted": sum(row["grade"] in ACCEPTED for row in nominated),
        "NominatedAssembled": sum(
            row["screen_overlap"] == "same_library_sylph"
            and int(row["assembly_support_bases"]) >= assembly_bases
            for row in overlap
        ),
        "AssemblySupportBases": assembly_bases,
        "OtherLibraryNominatedAccepted": sum(
            row["screen_overlap"] == "other_library_sylph" for row in overlap
        ),
        "ProkaryoticDetections": len(overlap),
        "NominatedBelow": sum(row["grade"] not in ACCEPTED for row in nominated),
        "NominatedTrace": sum(row["grade"] == "trace" for row in nominated),
        "NominatedNone": sum(row["grade"] == "none" for row in nominated),
        "MagOnlyGenomes": len(mag_only),
        "MagOnlyPresences": mag_presences,
        "MycoMagGenomes": len({row["target_id"] for row in myco}),
        "MycoMagDetections": len(myco),
        "MycoMagLibraries": len(myco_samples),
        "MycoBrackenLibraries": len(myco_samples & set(myco_labels["Bracken"])),
        "MycoPhyloflashLibraries": len(myco_samples & set(myco_labels["phyloFlash"])),
        "ArchaealMagGenomes": len({row["target_id"] for row in archaea}),
        "ArchaealMagDetections": len(archaea),
        "InsectsBefore": sum(n > 0 for n in before["Insecta"].values()),
        "InsectsAfter": sum(n > 0 for n in after["Insecta"].values()),
        "InsectMaxAfter": max(after["Insecta"].values()),
        "ReadThreshold": read_threshold,
        "Identity": identity,
        "EukDetections": len(euk_accepted),
        "EukLibraries": len({row["sample_id"] for row in euk_accepted}),
        "EukPhyloflashAssemblies": len(euk_phyloflash),
        "EukSharedAssemblies": len(euk_shared),
        "EukAssemblyOnly": len(euk_assembly_only),
        "KinetoplastidRescueReads": int(rescue[0]["read_pairs_ge97"]),
    }
    write("method_summary.tsv", method_rows)
    write("sylph_nomination_grades.tsv", nominated)
    write("ssu_control_comparison.tsv", control_rows)
    write("sylph_sources.tsv", nominations)
    write("genus_corroboration.tsv", corroboration)
    write("bacterial_detection_overlap.tsv", overlap)
    write("mycoplasmatales_screen_support.tsv", myco_support)
    write("eukaryote_detection_overlap.tsv", euk_accepted)
    (OUT / "summary.json").write_text(json.dumps(values, indent=2) + "\n")
    (OUT / "values.tex").write_text(
        "% Generated by scripts/summarize_screening_results.py; do not edit.\n"
        + "".join(
            "\\newcommand{\\screen" + key + "}{" + f"{value:,}" + "}\n"
            for key, value in values.items()
        )
    )
    snapshot = OUT / "provenance/summarize_screening_results.py"
    snapshot.parent.mkdir(exist_ok=True)
    shutil.copyfile(__file__, snapshot)
    outputs = [
        p for p in OUT.glob("*") if p.is_file() and p.name != "provenance.json"
    ] + [snapshot]
    if any(sha256(ROOT / name) != expected for name, expected in INPUTS.items()):
        raise ValueError("Source tables changed during reporting")
    provenance = {
        "status": "complete",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "command": "python3 scripts/summarize_screening_results.py",
        "source_base_git_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "script_sha256": sha256(Path(__file__)),
        "source_snapshot": str(snapshot.relative_to(ROOT)),
        "inputs_sha256": INPUTS,
        "outputs_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in outputs},
        "interpretation": "Taxonomic nominations are not validated presences. Coverage failures are not a measured classifier false-positive rate. SSU controls compare retained initial and host-aware assignments; no sequence analysis was rerun.",
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(values, indent=2))


if __name__ == "__main__":
    main()
