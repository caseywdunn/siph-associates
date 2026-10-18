#!/usr/bin/env python3
"""Audit COI reference origins and verify preserved candidates competitively.

Recruitment is deliberately unchanged. Suspect reference labels cannot establish
a crustacean assignment; divergent or conflicting matches remain candidates.
"""

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path: Path, rows: list[dict]) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def origin(lineage: list[str]) -> str:
    """Use NCBI lineage, allowing for names used before/after Crustacea changes."""
    crustaceans = {
        "Crustacea",
        "Pancrustacea",
        "Altocrustacea",
        "Oligostraca",
        "Malacostraca",
        "Copepoda",
        "Thecostraca",
        "Cirripedia",
        "Branchiopoda",
        "Remipedia",
        "Cephalocarida",
        "Ichthyostraca",
        "Ostracoda",
    }
    if "Hexapoda" not in lineage and crustaceans.intersection(lineage):
        return "crustacean"
    return "other_eukaryote"


def anchors(config: dict) -> None:
    root, source = Path(config["run"]), Path(config["source_panel"])
    out = root / "anchors"
    out.mkdir(parents=True, exist_ok=True)
    rows = table(source / "reference_metadata.tsv")
    annotated = {r["source_id"] for r in rows if r["source"] == "RefSeq"}
    lineages = {}
    with gzip.open(config["mitochondrial_genbank"], "rt") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            if record.id in annotated:
                lineages[record.id] = record.annotations.get("taxonomy", [])
    if annotated - lineages.keys():
        raise ValueError("Missing RefSeq lineage")
    for row in rows:
        row["origin"] = (
            origin(lineages[row["source_id"]])
            if row["source"] == "RefSeq"
            else "crustacean"
            if row["target"] == "1"
            else "other_eukaryote"
        )
        row["lineage"] = ";".join(lineages.get(row["source_id"], []))
    # Keep source target designations for provenance; origin corrects the scope.
    write_table(out / "reference_metadata.tsv", rows)
    metadata = {r["reference_id"]: r for r in rows}
    protein_rows = []
    with (out / "reference.faa").open("w") as handle:
        for record in SeqIO.parse(source / "annotated.faa", "fasta"):
            SeqIO.write(record, handle, "fasta")
            row = metadata[record.id]
            protein_rows.append(
                dict(
                    reference_id=record.id,
                    origin=row["origin"],
                    species=row["species"],
                    source_id=row["source_id"],
                )
            )
        for row in table(root / "source/uniprot_oxidases.tsv"):
            rid = "UP_" + row["Entry"]
            handle.write(f">{rid}\n{row['Sequence']}\n")
            protein_rows.append(
                dict(
                    reference_id=rid,
                    origin="prokaryote",
                    species=row["Organism"],
                    source_id=row["Entry"],
                )
            )
    write_table(out / "protein_metadata.tsv", protein_rows)
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                proteins_by_origin=Counter(r["origin"] for r in protein_rows),
                additional_crustacean_references=sum(
                    r["origin"] == "crustacean" and r["target"] == "0" for r in rows
                ),
            ),
            indent=2,
        )
        + "\n"
    )


def hits(path: Path) -> dict[str, list[dict]]:
    result = defaultdict(list)
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            # DIAMOND: qseqid sseqid pident length qcovhsp bitscore evalue
            result[row[0]].append(
                dict(
                    reference_id=row[1],
                    identity=float(row[2]),
                    length=int(row[3]),
                    coverage=float(row[4]),
                    bitscore=float(row[5]),
                )
            )
    return result


def deduplicate(config: dict) -> None:
    """Audit identical protein sequences once while retaining every source label."""
    root, source = Path(config["run"]), Path(config["source_panel"])
    representatives, mapping = {}, []
    (root / "anchors").mkdir(parents=True, exist_ok=True)
    with (root / "anchors/target_unique.faa").open("w") as handle:
        for record in SeqIO.parse(source / "target.faa", "fasta"):
            sequence = str(record.seq)
            if sequence not in representatives:
                representatives[sequence] = record.id
                SeqIO.write(record, handle, "fasta")
            mapping.append(
                dict(reference_id=record.id, representative=representatives[sequence])
            )
    write_table(root / "anchors/target_representatives.tsv", mapping)


def near_best(alignments: list[dict], fraction: float) -> list[dict]:
    if not alignments:
        return []
    best = max(h["bitscore"] for h in alignments)
    return [h for h in alignments if h["bitscore"] >= fraction * best]


def curate(config: dict) -> None:
    root, source = Path(config["run"]), Path(config["source_panel"])
    metadata = {
        r["reference_id"]: r for r in table(root / "anchors/protein_metadata.tsv")
    }
    manual = {r["reference_id"] for r in table(Path(config["manual_quarantine"]))}
    alignments = hits(root / "anchors/target_audit.tsv")
    representatives = {
        r["reference_id"]: r["representative"]
        for r in table(root / "anchors/target_representatives.tsv")
    }
    review, quarantine = [], set()
    for record in SeqIO.parse(source / "target.faa", "fasta"):
        near = near_best(
            alignments.get(representatives[record.id], []),
            config["reference_near_best_fraction"],
        )
        origins = {metadata[h["reference_id"]]["origin"] for h in near}
        reason = ""
        if record.id in manual:
            reason = "independent_source_record_and_bacterial_homology_audit"
        elif "prokaryote" in origins:
            reason = "prokaryotic_oxidase_among_best_anchors"
        elif origins != {"crustacean"}:
            reason = "crustacean_origin_not_resolved_by_annotated_anchors"
        if reason:
            quarantine.add(record.id)
        review.append(
            dict(
                reference_id=record.id,
                assignment_status="quarantine" if reason else "retained",
                reason=reason,
                anchor_origins=";".join(sorted(origins)),
                best_bitscore=max((h["bitscore"] for h in near), default=0),
                near_best_anchors=";".join(sorted({h["reference_id"] for h in near})),
            )
        )
    out = root / "curated"
    out.mkdir(parents=True, exist_ok=True)
    write_table(out / "reference_audit.tsv", review)
    source_rows = table(root / "anchors/reference_metadata.tsv")
    for row in source_rows:
        row["assignment_status"] = (
            "quarantine" if row["reference_id"] in quarantine else "retained"
        )
        if row["reference_id"] in quarantine:
            row["origin"] = "unresolved"
    write_table(out / "reference_metadata.tsv", source_rows)
    held_ids = {
        r["reference_id"]
        for r in table(Path(config["validation_run"]) / "panel/reference_metadata.tsv")
    }
    used, protein_rows = set(), []
    with (out / "full.faa").open("w") as full, (out / "heldout.faa").open("w") as held:
        for path in [root / "anchors/reference.faa", source / "target.faa"]:
            for record in SeqIO.parse(path, "fasta"):
                if record.id in used or record.id in quarantine:
                    continue
                used.add(record.id)
                row = metadata.get(record.id)
                if row is None:
                    row = dict(
                        reference_id=record.id,
                        origin="crustacean",
                        species="",
                        source_id=record.id,
                    )
                protein_rows.append(row)
                SeqIO.write(record, full, "fasta")
                if record.id.startswith("UP_") or record.id in held_ids:
                    SeqIO.write(record, held, "fasta")
    write_table(out / "protein_metadata.tsv", protein_rows)
    SeqIO.write(
        (r for r in SeqIO.parse(source / "target.fna", "fasta") if r.id in quarantine),
        out / "quarantined.fna",
        "fasta",
    )
    (out / "summary.json").write_text(
        json.dumps(
            dict(
                audited=len(review),
                quarantined=len(quarantine),
                quarantine_reasons=Counter(r["reason"] for r in review if r["reason"]),
                manual_quarantine_included=sorted(manual & quarantine),
                assignment_proteins_by_origin=Counter(
                    r["origin"] for r in protein_rows
                ),
                recruitment="Original broad panel unchanged; all recruited pairs retained",
            ),
            indent=2,
        )
        + "\n"
    )


def queries(config: dict) -> None:
    root, pilot, validation = map(
        Path, (config["run"], config["pilot_run"], config["validation_run"])
    )
    out = root / "queries"
    out.mkdir(parents=True, exist_ok=True)
    for dataset in ["pilot", "controls", "heldout"]:
        directories = (
            sorted(pilot.glob("Ahuja2026__*"))
            if dataset == "pilot"
            else [
                pilot / "controls" if dataset == "controls" else validation / "heldout"
            ]
        )
        pairs = []
        with (
            (out / f"{dataset}.fna").open("w") as fasta,
            (out / f"{dataset}.nt.tsv").open("w") as nt,
        ):
            for directory in directories:
                prefix = directory.name + "__"
                for record in SeqIO.parse(directory / "candidates.fna", "fasta"):
                    fasta.write(f">{prefix}{record.id}\n{record.seq}\n")
                for row in table(directory / "candidate_pairs.tsv"):
                    pairs.append(
                        dict(
                            dataset=dataset,
                            sample=directory.name,
                            **row,
                            query_prefix=prefix + row["candidate_id"],
                        )
                    )
                with (directory / "competitive.tsv").open() as handle:
                    for line in handle:
                        nt.write(prefix + line)
        write_table(out / f"{dataset}.pairs.tsv", pairs)


def protein_evidence(
    alignments: list[dict], metadata: dict, config: dict
) -> tuple[str, set[str]]:
    near = near_best(alignments, config["read_near_best_fraction"])
    origins = {metadata[h["reference_id"]]["origin"] for h in near}
    if not near:
        return "unresolved", set()
    if len(origins) > 1:
        return "ambiguous", origins
    if any(
        h["identity"] < config["minimum_protein_identity"]
        or h["length"] < config["minimum_protein_alignment_aa"]
        or h["coverage"] < config["minimum_protein_query_coverage"]
        for h in near
    ):
        return "weak", origins
    return next(iter(origins)), origins


def nucleotide_evidence(alignments: list[dict], metadata: dict, config: dict) -> str:
    near = near_best(alignments, config["read_near_best_fraction"])
    if not near:
        return "unresolved"
    origins = {metadata[h["reference_id"]]["origin"] for h in near}
    if len(alignments) >= 200 and min(h["bitscore"] for h in alignments) >= config[
        "read_near_best_fraction"
    ] * max(h["bitscore"] for h in alignments):
        return "ambiguous"
    if "unresolved" in origins or len(origins) > 1:
        return "ambiguous"
    if any(
        h["identity"] < config["minimum_nucleotide_identity"]
        or h["length"] < config["minimum_nucleotide_alignment_nt"]
        for h in near
    ):
        return "weak"
    return next(iter(origins))


def decide_pair(proteins: list[str], nucleotides: list[str]) -> str:
    """Require agreement for eukaryotic calls; keep weaker evidence for review."""
    informative = set(proteins) & {"crustacean", "other_eukaryote", "prokaryote"}
    if "ambiguous" in proteins or len(informative) > 1:
        return "ambiguous"
    if informative == {"prokaryote"}:
        if set(nucleotides) & {"crustacean", "other_eukaryote"}:
            return "ambiguous"
        return "bacterial_or_archaeal_oxidase"
    for label, result in [
        ("crustacean", "crustacean_supported"),
        ("other_eukaryote", "other_eukaryote_supported"),
    ]:
        if informative == {label} and label in nucleotides:
            opposite = "other_eukaryote" if label == "crustacean" else "crustacean"
            if opposite in nucleotides or "ambiguous" in nucleotides:
                return "ambiguous"
            # Support must be present in the same mate, not unrelated mate matches.
            if any(p == n == label for p, n in zip(proteins, nucleotides)):
                return result
    return "unresolved_candidate"


def report(config: dict) -> None:
    root = Path(config["run"])
    protein_meta = {
        r["reference_id"]: r for r in table(root / "curated/protein_metadata.tsv")
    }
    nt_meta = {
        r["reference_id"]: r for r in table(root / "curated/reference_metadata.tsv")
    }
    all_rows, summaries, gates = [], {}, {}
    for dataset in ["pilot", "controls", "heldout"]:
        protein = defaultdict(list)
        for code in [1, 5]:
            for query, rows in hits(
                root / "verification" / dataset / f"code{code}.tsv"
            ).items():
                protein[query].extend(rows)
        nt = defaultdict(list)
        with (root / "queries" / f"{dataset}.nt.tsv").open() as handle:
            for r in csv.reader(handle, delimiter="\t"):
                nt[r[0]].append(
                    dict(
                        reference_id=r[1],
                        identity=float(r[2]),
                        length=int(r[3]),
                        bitscore=float(r[8]),
                    )
                )
        details = []
        for row in table(root / "queries" / f"{dataset}.pairs.tsv"):
            queries = [row["query_prefix"] + f"_{mate}" for mate in [1, 2]]
            aa = [protein_evidence(protein[q], protein_meta, config) for q in queries]
            dna = [nucleotide_evidence(nt[q], nt_meta, config) for q in queries]
            best_refs = sorted(
                {
                    h["reference_id"]
                    for q in queries
                    for h in near_best(protein[q], config["read_near_best_fraction"])
                }
            )
            details.append(
                dict(
                    **row,
                    classification=decide_pair([a[0] for a in aa], dna),
                    protein_R1=aa[0][0],
                    protein_R2=aa[1][0],
                    nucleotide_R1=dna[0],
                    nucleotide_R2=dna[1],
                    protein_origins=";".join(sorted(aa[0][1] | aa[1][1])),
                    near_best_proteins=";".join(best_refs),
                )
            )
        summaries[dataset] = dict(
            pairs=len(details),
            classifications=Counter(r["classification"] for r in details),
        )
        if dataset != "pilot":
            truth_path = (
                Path(config["pilot_run"]) / "controls/truth.json"
                if dataset == "controls"
                else Path(config["validation_run"]) / "heldout/truth.json"
            )
            truth = json.loads(truth_path.read_text())
            found = {r["read_pair_id"]: r["classification"] for r in details}
            positive = [r for r in truth if r["expected"] == "crustacean"]
            negative = [r for r in truth if r["expected"] != "crustacean"]
            summaries[dataset]["positive_classifications"] = Counter(
                found.get(r["read_pair_id"], "not_recruited") for r in positive
            )
            summaries[dataset]["negative_classifications"] = Counter(
                found.get(r["read_pair_id"], "not_recruited") for r in negative
            )
            if dataset == "controls":
                gates["original_positive_controls_include_supported_match"] = any(
                    found.get(r["read_pair_id"]) == "crustacean_supported"
                    for r in positive
                )
            gates[f"{dataset}_positives_recruited_and_not_rejected"] = all(
                found.get(r["read_pair_id"])
                in {"crustacean_supported", "ambiguous", "unresolved_candidate"}
                for r in positive
            )
            gates[f"{dataset}_negatives_not_called_crustacean"] = all(
                found.get(r["read_pair_id"]) != "crustacean_supported" for r in negative
            )
            by_group = defaultdict(Counter)
            for row in truth:
                by_group[row.get("group", row["expected"])][
                    found.get(row["read_pair_id"], "not_recruited")
                ] += 1
            summaries[dataset]["by_group"] = by_group
        else:
            known_bacterial = [
                r
                for r in details
                if r["sample"] == "Ahuja2026__NA35"
                and r["candidate_id"] in {"c00000008", "c00000016"}
            ]
            gates["both_known_bacterial_pairs_recognized"] = len(
                known_bacterial
            ) == 2 and all(
                r["classification"] == "bacterial_or_archaeal_oxidase"
                for r in known_bacterial
            )
            gates["all_92_pilot_pairs_preserved"] = len(details) == 92
        all_rows.extend(details)
    curation = json.loads((root / "curated/summary.json").read_text())
    gates["four_flagged_barcodes_quarantined"] = (
        len(curation["manual_quarantine_included"]) == 4
    )
    summaries["reference_audit"] = curation
    summaries["gates"] = gates
    summaries["passed"] = all(gates.values())
    summaries["limitations"] = (
        "Control retention tests reference divergence, not all natural prey or all sequencing errors. Broad recruitment retains false candidates. Short conserved COI often remains unresolved. Supported calls are marker evidence, not species or dietary assignments."
    )
    (root / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    write_table(root / "classified_pairs.tsv", all_rows)
    (root / "gate.json").write_text(
        json.dumps(dict(passed=all(gates.values()), checks=gates), indent=2) + "\n"
    )
    paths = [
        root / "summary.json",
        root / "classified_pairs.tsv",
        root / "gate.json",
        root / "curated/reference_audit.tsv",
        root / "curated/reference_metadata.tsv",
        root / "curated/protein_metadata.tsv",
    ]
    (root / "output_sha256.txt").write_text(
        "".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p}\n" for p in paths)
    )
    print(json.dumps(dict(summaries=summaries, passed=all(gates.values())), indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage", choices=["anchors", "deduplicate", "curate", "queries", "report"]
    )
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    globals()[args.stage](json.loads(args.config.read_text()))


if __name__ == "__main__":
    main()
