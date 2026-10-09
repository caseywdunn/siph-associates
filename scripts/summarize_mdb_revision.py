#!/usr/bin/env python3
"""Summarize saved detections and annotations for the MDB manuscript revision.

Standard library only. Does not rerun annotation, mapping, or hypothesis tests.
See docs/mdb_revision.md for scope, denominators, and interpretation limits.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
import subprocess
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from statistics import median

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/mdb_revision"
GENES = "data/results/mycoplasmatales_gene_content/summary"
GENERA = ("Vibrio", "Pseudoalteromonas", "Alteromonas", "Photobacterium")
FOCAL = ("MAGSP0005", "MAGSP0007", "MAGSP0010", "MAGSP0011", "MAGSP0012")
GROUPS = {
    "Metamycoplasmataceae": "metamycoplasmataceae_detected",
    "clade A": "clade_a_detected",
    "clade B": "clade_b_detected",
    "DT-68": "dt68_detected",
}
INPUTS: dict[str, str] = {}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(relative: str) -> list[dict[str, str]]:
    path = ROOT / relative
    INPUTS[relative] = sha256(path)
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name: str, rows: list[dict]) -> None:
    with (OUT / name).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict], scope: str, stratum: str = "") -> list[dict]:
    results = []
    for group, column in GROUPS.items():
        positive = [r for r in rows if int(r[column])]
        negative = [r for r in rows if not int(r[column])]
        for partner in (*GENERA, "any_candidate_genus"):
            both = sum(int(r[partner]) for r in positive)
            partner_only = sum(int(r[partner]) for r in negative)
            results.append(
                {
                    "scope": scope,
                    "stratum": stratum,
                    "group": group,
                    "partner": partner,
                    "n": len(rows),
                    "group_positive": len(positive),
                    "group_negative": len(negative),
                    "both": both,
                    "group_only": len(positive) - both,
                    "partner_only": partner_only,
                    "neither": len(negative) - partner_only,
                    "partner_fraction_group_positive": both / len(positive)
                    if positive
                    else "",
                    "partner_fraction_group_negative": partner_only / len(negative)
                    if negative
                    else "",
                    "median_read_pairs_group_positive": median(
                        int(r["read_pairs"]) for r in positive
                    )
                    if positive
                    else "",
                    "median_read_pairs_group_negative": median(
                        int(r["read_pairs"]) for r in negative
                    )
                    if negative
                    else "",
                }
            )
    return results


def cooccurrence() -> None:
    context = read("figures/collection_context/library_context.tsv")
    samples = {r["sample_id"]: dict(r) for r in context}
    if len(samples) != len(context) or len({r["library_id"] for r in context}) != len(
        context
    ):
        raise ValueError("Duplicate sample or library identifier")
    detections = read("figures/cohort/bacterial_detections.tsv")
    seen = set()
    family_samples = set()
    for r in samples.values():
        r.update(dict.fromkeys(GENERA, 0))
    for row in detections:
        key = (row["sample_id"], row["target_id"])
        if key in seen or row["sample_id"] not in samples:
            raise ValueError(f"Duplicate or unjoined detection: {key}")
        seen.add(key)
        if row["grade"] not in {"validated", "high_confidence"}:
            raise ValueError("The export must contain accepted detections only")
        taxonomy = row["taxonomy"].split(";")
        if "f__Metamycoplasmataceae" in taxonomy:
            family_samples.add(row["sample_id"])
        for genus in GENERA:
            if f"g__{genus}" in taxonomy:
                samples[row["sample_id"]][genus] = 1
    if family_samples != {
        r["sample_id"] for r in context if int(r[GROUPS["Metamycoplasmataceae"]])
    }:
        raise ValueError("Family detections disagree with collection-context summary")
    for r in samples.values():
        r["any_candidate_genus"] = int(any(r[g] for g in GENERA))
        if r["Vibrio"] != int(r["vibrio_detected"]):
            raise ValueError("Vibrio detections disagree with existing summary")
    rows = list(samples.values())
    write("specimens.tsv", rows)
    physalia = [r for r in rows if r["host_genus"] == "Physalia"]
    eligible = [r for r in physalia if r["single_flowcell_eligible"] == "true"]
    scopes = {
        "all_specimens": rows,
        "Physalia": physalia,
        "Physalia_single_flowcell": eligible,
        "Physalia_tentacle": [r for r in physalia if r["tissue_group"] == "tentacle"],
        "P_utriculus_tentacle": [
            r
            for r in physalia
            if r["tissue_group"] == "tentacle"
            and r["species_current"] == "Physalia utriculus"
        ],
    }
    summary = [
        row for scope, subset in scopes.items() for row in summarize(subset, scope)
    ]
    write("cooccurrence.tsv", summary)
    strata = defaultdict(list)
    for r in eligible:
        strata[
            tuple(
                r[k] or "unknown"
                for k in (
                    "species_current",
                    "ocean_region",
                    "tissue_group",
                    "flowcells",
                )
            )
        ].append(r)
    stratified = [
        row
        for key, subset in sorted(strata.items())
        for row in summarize(subset, "Physalia_single_flowcell", " | ".join(key))
    ]
    write("cooccurrence_strata.tsv", stratified)
    table = [
        r"\begin{table}[p]\centering\small",
        r"\begin{tabular}{lrr}\toprule",
        r"Candidate genus & With Metamycoplasmataceae & Without \\\midrule",
    ]
    for row in summary:
        if row["scope"] != "Physalia" or row["group"] != "Metamycoplasmataceae":
            continue
        label = (
            "Any of these genera"
            if row["partner"] == "any_candidate_genus"
            else r"\textit{" + row["partner"] + "}"
        )
        table.append(
            f"{label} & {row['both']}/{row['group_positive']} & {row['partner_only']}/{row['group_negative']} "
            + r"\\"
        )
    table += [
        r"\bottomrule\end{tabular}",
        r"\caption{Descriptive co-detection in \textit{Physalia}. Fractions give specimens with at least one accepted detection in each candidate genus over specimens with or without Metamycoplasmataceae. Each specimen is counted once per row; genera are not mutually exclusive. Genus labels do not establish chitin degradation in the detected strains. These unadjusted frequencies do not test cross-feeding or control geographic, host, tissue, sequencing-depth or batch differences. Supporting tables retain clade-specific counts, restricted subsets and strata defined by host species, region, tissue and physical flowcell.}",
        r"\label{tab:cooccurrence}\end{table}",
    ]
    (OUT / "cooccurrence_table.tex").write_text("\n".join(table) + "\n")
    combined = next(
        r
        for r in summary
        if r["scope"] == "Physalia"
        and r["group"] == "Metamycoplasmataceae"
        and r["partner"] == "any_candidate_genus"
    )
    (OUT / "cooccurrence_results.tex").write_text(
        "In \\textit{Physalia}, at least one of the candidate genera occurred in\n"
        f"{combined['both']} of {combined['group_positive']} specimens with Metamycoplasmataceae and\n"
        f"{combined['partner_only']} of {combined['group_negative']} specimens without them\n"
        "(Table~\\ref{tab:cooccurrence}). The genus-specific unadjusted frequencies\n"
        "varied in direction and did not show a consistent increase with\n"
        "Metamycoplasmataceae. These descriptive comparisons do not establish\n"
        "nutritional exchange or an adjusted ecological association.\n"
    )


def gene_content() -> None:
    stats = read(f"{GENES}/genome_stats.tsv")
    quality = {r["genome"]: r for r in stats if r["genome"].startswith("MAGSP")}
    if {g for g, r in quality.items() if r["near_complete"].lower() == "true"} != set(
        FOCAL
    ):
        raise ValueError(
            "Near-complete focal genome set changed; review reporting scope"
        )
    focal = read(f"{GENES}/focal_functions.tsv")
    modules = {
        tier: read(f"{GENES}/module_completeness{suffix}.tsv")
        for tier, suffix in (("strict", "_strict"), ("relaxed", ""))
    }
    # These are the stored biosynthetic modules, not a complete auxotrophy screen.
    names = {
        r["module"]: r["name"]
        for rs in modules.values()
        for r in rs
        if "Amino acid" in r["class"]
        and "biosynthesis" in r["name"]
        and "plants" not in r["name"]
        and "Ethylene" not in r["name"]
    }
    records = []
    lookup = {}
    for tier, rs in modules.items():
        by_pair = {(r["genome"], r["module"]): r for r in rs}
        for genome, q in quality.items():
            for module, name in sorted(names.items()):
                row = by_pair.get((genome, module))
                value = row["completeness"] if row else ""
                lookup[(tier, genome, module)] = value
                records.append(
                    {
                        "genome": genome,
                        "lineage": q["lineage"],
                        "near_complete": q["near_complete"],
                        "tier": tier,
                        "module": module,
                        "name": name,
                        "completeness": value,
                        "status": "reported" if row else "not_returned",
                        "missing_ko": row["missing_ko"] if row else "",
                    }
                )
    write("amino_acid_modules.tsv", records)
    selected = [
        r
        for r in focal
        if any(
            term in r["function"]
            for term in ("CRISPR", "restriction-modification", "chitinase")
        )
    ]
    write(
        "defense_chitin_genes.tsv",
        [
            {
                k: v
                for k, v in r.items()
                if k in ("function", "ko", "definition", "assignment") or k in quality
            }
            for r in selected
        ],
    )
    write("genome_quality.tsv", list(quality.values()))
    genes = {(r["assignment"], r["ko"]): r for r in selected}
    table = [
        r"\begin{table}[p]\centering\footnotesize",
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{lrrrrr}\toprule",
        r"Feature & 0005 (A) & 0007 (B) & 0010 (B) & 0011 (DT-68) & 0012 (K) \\\midrule",
        r"\multicolumn{6}{l}{Selected amino-acid modules: completeness (\%)} \\",
    ]
    for module, label in (
        ("M00609", "Cysteine from methionine"),
        ("M00017", "Methionine from aspartate"),
        ("M00026", "Histidine from PRPP"),
        ("M00844", "Arginine from ornithine"),
    ):
        cells = []
        for genome in FOCAL:
            values = [lookup[(t, genome, module)] for t in ("strict", "relaxed")]
            cells.append("/".join(f"{float(v):.1f}" if v else "NR" for v in values))
        table.append(f"{label} & " + " & ".join(cells) + r" \\")
    table.append(
        r"\midrule\multicolumn{6}{l}{Selected defense genes: assigned copies} \\"
    )
    for ko, label in (
        ("K15342", "Cas1"),
        ("K09951", "Cas2"),
        ("K09952", "Cas9/Csn1"),
        ("K19137", "Csn2"),
        ("K01153", "Type I restriction R"),
        ("K03427", "Type I modification M"),
        ("K01154", "Type I specificity S"),
        ("K07316", "DNA methyltransferase"),
        ("K01156", "Type III restriction"),
    ):
        cells = [
            "/".join(genes[(t, ko)][g] for t in ("strict", "relaxed")) for g in FOCAL
        ]
        table.append(f"{label} & " + " & ".join(cells) + r" \\")
    table += [
        r"\bottomrule\end{tabular}",
        r"\caption{Selected nutritional and defense annotations in near-complete Metamycoplasmataceae genomes. Column identifiers omit the MAGSP prefix; A, B, DT-68 and K identify the four focal lineages. Cells give strict/relaxed assignments. NR means no row was returned in the stored module output and is not an imputed completeness estimate. Module identifiers are M00609, M00017, M00026 and M00844, respectively; these overlapping and precursor-dependent routes cannot be counted as experimentally established auxotrophies. Defense entries are selected orthologs, not counts of functional systems. Incomplete assemblies can lose genes. Machine-readable tables additionally retain the two presence-only genomes, all selected amino-acid module records, missing orthologs, and the screened chitin-degrading genes.}",
        r"\label{tab:nutrition-defense}\end{table}",
    ]
    (OUT / "nutrition_defense_table.tex").write_text("\n".join(table) + "\n")


def phage_summary() -> None:
    rows = read("data/results/phase6_analysis/v2/primary/phage_host_links.tsv")
    if len({r["votu_id"] for r in rows}) != len(rows):
        raise ValueError("Multiple host predictions per vOTU need revised wording")
    photo = [r for r in rows if "s__Photobacterium angustum" in r["host_taxonomy"]]
    others = [r for r in rows if r not in photo]
    if any(int(r["co_present_libraries"]) for r in photo) or not all(
        int(r["co_present_libraries"]) for r in others
    ):
        raise ValueError("Phage co-detection pattern changed; review interpretation")
    words = {
        0: "zero",
        1: "one",
        2: "two",
        3: "three",
        4: "four",
        5: "five",
        6: "six",
        7: "seven",
    }

    def n(value: int) -> str:
        return words.get(value, str(value))

    text = (
        "Bacterial CRISPR systems retain fragments of viral DNA from previous\n"
        f"encounters. {n(sum(int(r['spacers']) for r in rows)).capitalize()} matches to these fragments supported host predictions\n"
        f"for {n(len(rows))} viral groups. {n(len(photo)).capitalize()} had \\textit{{Photobacterium angustum}} as their\n"
        "predicted host, but none of these viral groups was detected in the same\n"
        "specimen as that bacterium. The other " + n(len(others)) + " matched other\n"
        "gammaproteobacteria and were detected in the same specimens as their\n"
        "respective predicted hosts. These matches link part of the viral diversity\n"
        "to the associated bacterial community; spacer matches record earlier\n"
        "exposure and do not by themselves demonstrate infection in the sampled specimen.\n"
    )
    (OUT / "phage_summary.tex").write_text(text)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cooccurrence()
    gene_content()
    phage_summary()
    source = Path(__file__).resolve()
    manifest = {
        "status": "complete",
        "created_utc": datetime.now(UTC).isoformat(),
        "command": "python3 scripts/summarize_mdb_revision.py",
        "python": platform.python_version(),
        "source_base_git_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "script_sha256": sha256(source),
        "inputs_sha256": INPUTS,
        "outputs_sha256": {
            str(p.relative_to(ROOT)): sha256(p)
            for p in sorted(OUT.iterdir())
            if p.suffix in {".tsv", ".tex"}
        },
        "scope": "Descriptive summaries of saved results; no new significance tests, reannotation or sequence processing.",
    }
    (OUT / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "status": "complete",
                "output": str(OUT),
                "files": len(manifest["outputs_sha256"]),
            }
        )
    )


if __name__ == "__main__":
    main()
