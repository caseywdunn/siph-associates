#!/usr/bin/env python3
"""Render publication figures from retained Mycoplasmatales results.

Run from the repository root on a compute node. This does not infer a new tree,
reassign functions, or estimate missing genes. Source tables and input checksums
are retained beside the figures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import Bio
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

FUNCTIONS = [
    ("K01478", "Arginine deiminase (arcA)"),
    ("K00611", "Ornithine carbamoyltransferase (arcB)"),
    ("K00926", "Carbamate kinase (arcC)"),
    ("K00016", "L-lactate dehydrogenase"),
    ("K00166", "2-oxoacid dehydrogenase E1 alpha"),
    ("K00167", "2-oxoacid dehydrogenase E1 beta"),
    ("K00627", "Pyruvate dehydrogenase E2"),
    ("K00382", "Dihydrolipoyl dehydrogenase E3"),
    ("K00625", "Phosphotransacetylase (pta)"),
    ("K00925", "Acetate kinase (ackA)"),
    ("K02440", "Glycerol uptake facilitator (glpF)"),
    ("K00864", "Glycerol kinase (glpK)"),
    ("K00105", "Glycerophosphate oxidase (glpO)"),
    ("K01639", "N-acetylneuraminate lyase (nanA)"),
    ("K01788", "N-acylglucosamine-P epimerase (nanE)"),
    ("K01186", "Sialidase"),
]
MAG_ORDER = [
    "MAGSP0005", "MAGSP0007", "MAGSP0010", "MAGSP0029",
    "MAGSP0031", "MAGSP0011", "MAGSP0012",
]
COLORS = {"not_detected": "#EDF0F2", "relaxed": "#EBA34E", "strict": "#197F8C"}


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a retained tabular result without changing its content."""
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_rows(path: Path, rows: list[dict]) -> None:
    """Write the plotted values so readers can inspect every mark."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def checksum(path: Path) -> str:
    """Calculate a SHA-256 identity for a small input or figure."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_figure(fig: plt.Figure, output: Path, name: str) -> None:
    """Save vector PDF with embedded TrueType fonts and a review preview."""
    fig.savefig(output / f"{name}.pdf", facecolor="white", bbox_inches="tight")
    fig.savefig(output / f"{name}.png", dpi=160, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def gene_functions(focal: list[dict], genomes: list[dict], output: Path) -> None:
    """Plot assignments rather than inferring phenotype or gene absence."""
    lookup = {(row["ko"], row["assignment"]): row for row in focal}
    if len(lookup) != len(focal):
        raise ValueError("Duplicate KO/tier records in focal_functions.tsv")
    quality = {row["genome"]: row for row in genomes}
    family = [row["genome"] for row in genomes if row["genome"].startswith("GTDB_")
              and row["near_complete"] == "True" and row["included"] == "True"]
    if len(family) != 130 or not all(mag in quality for mag in MAG_ORDER):
        raise ValueError("Genome set changed: inspect the selected figure scope")
    matrix = np.zeros((len(FUNCTIONS), len(MAG_ORDER)), dtype=int)
    family_counts = {tier: [] for tier in ("strict", "relaxed")}
    plotted = []
    for i, (ko, label) in enumerate(FUNCTIONS):
        strict, relaxed = lookup[ko, "strict"], lookup[ko, "relaxed"]
        for tier in family_counts:
            family_counts[tier].append(sum(int(lookup[ko, tier][g]) > 0 for g in family))
        for j, genome in enumerate(MAG_ORDER):
            counts = {"strict": int(strict[genome]), "relaxed": int(relaxed[genome])}
            if counts["strict"] > counts["relaxed"]:
                raise ValueError(f"Strict counts exceed relaxed counts: {ko}, {genome}")
            matrix[i, j] = 2 if counts["strict"] else int(counts["relaxed"] > 0)
            plotted.append({"ko": ko, "label": label, "genome": genome,
                            "strict_copies": counts["strict"], "relaxed_copies": counts["relaxed"],
                            "near_complete": quality[genome]["near_complete"],
                            "family_n": len(family),
                            "family_strict_present": family_counts["strict"][-1],
                            "family_relaxed_present": family_counts["relaxed"][-1]})
    write_rows(output / "source_gene_functions.tsv", plotted)

    fig = plt.figure(figsize=(8.0, 6.35))
    grid = fig.add_gridspec(1, 2, width_ratios=[7, 2.3], left=0.365,
                           right=0.965, bottom=0.22, top=0.82, wspace=0.2)
    ax, frequency = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    ax.imshow(matrix, cmap=ListedColormap(list(COLORS.values())), vmin=0, vmax=2,
              aspect="auto", interpolation="none")
    ax.set_yticks(range(len(FUNCTIONS)), [label for _, label in FUNCTIONS], fontsize=8)
    ax.set_xticks(range(len(MAG_ORDER)),
                  [g + ("*" if quality[g]["near_complete"] != "True" else "")
                   for g in MAG_ORDER], rotation=60, ha="right", fontsize=8)
    ax.tick_params(length=0)
    ax.set_xticks(np.arange(-0.5, len(MAG_ORDER), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(FUNCTIONS), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.2)
    ax.tick_params(which="minor", length=0)
    for j, text in [(0, "Clade A"), (2.5, "Clade B"), (5, "DT-68"), (6, "Myco. K")]:
        ax.text(j, -1.12, text, ha="center", va="center", fontsize=8, clip_on=False)
    for left, right in [(-0.4, 0.4), (0.6, 4.4), (4.6, 5.4), (5.6, 6.4)]:
        ax.plot([left, right], [-0.78, -0.78], color="#68777B", lw=1, clip_on=False)
    for i in range(len(FUNCTIONS)):
        for j, genome in enumerate(MAG_ORDER):
            copies = int(lookup[FUNCTIONS[i][0], "relaxed"][genome])
            if copies > 1:
                ax.text(j, i, str(copies), ha="center", va="center", fontsize=8,
                        color="white" if matrix[i, j] == 2 else "#333333")
    y = np.arange(len(FUNCTIONS))
    frequency.barh(y - 0.16, np.asarray(family_counts["strict"]) * 100 / len(family),
                   height=0.29, color=COLORS["strict"])
    frequency.barh(y + 0.16, np.asarray(family_counts["relaxed"]) * 100 / len(family),
                   height=0.29, color=COLORS["relaxed"])
    frequency.set(ylim=ax.get_ylim(), xlim=(0, 100), yticks=[], xticks=[0, 50, 100])
    frequency.set_xlabel("Genomes (%)", fontsize=8)
    frequency.set_title("Family context\n(n = 130)", fontsize=9, pad=10)
    frequency.tick_params(labelsize=8, length=3)
    frequency.spines[["top", "right", "left"]].set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for boundary in [2.5, 3.5, 9.5, 12.5]:
        ax.axhline(boundary, color="#87969A", lw=0.8)
    fig.text(0.025, 0.965, "Predicted functions in siphonophore Mycoplasmatales",
             fontsize=12, weight="bold", va="top")
    fig.text(0.025, 0.915, "Conserved ortholog assignments; strict and relaxed thresholds",
             fontsize=9, color="#48565C")
    legend = [Patch(facecolor=COLORS["strict"], label="Strict (also relaxed)"),
              Patch(facecolor=COLORS["relaxed"], label="Relaxed only"),
              Patch(facecolor=COLORS["not_detected"], label="Not detected")]
    fig.legend(handles=legend, loc="lower left", bbox_to_anchor=(0.02, 0.035),
               ncol=3, frameon=False, fontsize=8)
    fig.text(0.025, 0.025,
             "* Presence-only MAGs: incomplete or contamination-flagged. Numbers: relaxed-tier copies >1.\n"
             "Family bars show strict (teal) and all relaxed (orange) assignments; no phenotype is measured.",
             fontsize=7.5, va="top")
    save_figure(fig, output, "gene_functions")


def node_support(clade) -> float | None:
    """Read FastTree support, including GTDB's support:taxonomy node labels."""
    if clade.confidence is not None:
        return float(clade.confidence)
    if clade.name:
        try:
            return float(clade.name.split(":", 1)[0])
        except ValueError:
            return None
    return None


def focus_tree(
    tree_path: Path, genomes: list[dict], labels: list[dict], output: Path
) -> None:
    """Display the existing MRCA subtree without removing duplicate source tips."""
    tree = Phylo.read(tree_path, "newick")
    terminals = tree.get_terminals()
    names = [tip.name for tip in terminals]
    if len(names) != len(set(names)) or not set(MAG_ORDER).issubset(names):
        raise ValueError("Tree tips are duplicated or missing selected study MAGs")
    root = tree.common_ancestor(MAG_ORDER)
    tips = root.get_terminals()
    if len(terminals) != 483 or len(tips) != 21:
        raise ValueError(
            "Tree changed: inspect the selected 21-tip clade before plotting"
        )
    taxonomy = {row["genome"]: row for row in genomes}
    external = {row["genome"]: row for row in labels}
    depth, height = {}, {}

    def layout(clade, distance=0.0):
        depth[clade] = distance
        for child in clade.clades:
            layout(child, distance + (child.branch_length or 0.0))
        height[clade] = (
            (len(tips) - 1 - tips.index(clade))
            if clade.is_terminal()
            else (height[clade.clades[0]] + height[clade.clades[-1]]) / 2
        )

    layout(root)
    max_depth = max(depth.values())
    fig, ax = plt.subplots(figsize=(8.0, 7.8))
    fig.subplots_adjust(left=0.035, right=0.97, bottom=0.18, top=0.84)
    rows = []
    source_colors = {"study": "#197F8C", "external": "#B75B24", "GTDB": "#48565C"}
    mag_lineages = {
        "MAGSP0005": "clade A",
        "MAGSP0007": "clade B",
        "MAGSP0010": "clade B",
        "MAGSP0029": "clade B",
        "MAGSP0031": "clade B",
        "MAGSP0011": "DT-68",
        "MAGSP0012": "Mycoplasma_K",
    }
    for clade in root.find_clades(order="preorder"):
        if clade.clades:
            ax.plot(
                [depth[clade]] * 2,
                [height[clade.clades[0]], height[clade.clades[-1]]],
                color="#65767B",
                linewidth=0.8,
            )
            for child in clade.clades:
                ax.plot(
                    [depth[clade], depth[child]],
                    [height[child]] * 2,
                    color="#65767B",
                    linewidth=0.8,
                )
            support = node_support(clade)
            if support is not None:
                ax.text(
                    depth[clade] - 0.008 * max_depth,
                    height[clade] + 0.16,
                    f"{support:.3g}",
                    fontsize=9,
                    ha="right",
                    va="bottom",
                    color="#A65120" if support < 0.9 else "#37464A",
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.35},
                )
            continue
        tip = clade.name
        duplicate = "GCA_929200685.1" in tip
        if tip.startswith("MAGSP"):
            source = "study"
            label = f"{tip}  |  {mag_lineages[tip]}"
            classification = taxonomy[tip]["classification"]
        elif tip.startswith("EXT_"):
            source = "external"
            host = external[tip]["host_or_source"].split(";")[0].split(" (")[0]
            label = f"{host}  |  {tip.removeprefix('EXT_')}"
            classification = external[tip]["classification"]
        else:
            source = "GTDB"
            accession = tip.removeprefix("RS_").removeprefix("GB_")
            row = taxonomy.get(f"GTDB_{accession}") or taxonomy.get(f"EXT_{accession}")
            if row is None:
                raise ValueError(f"Missing retained taxonomy for {tip}")
            classification = row["classification"]
            label = classification.split(";")[-1].removeprefix("s__") or accession
        if duplicate:
            label += " [dup.]"
        color = source_colors[source]
        ax.scatter(
            depth[clade],
            height[clade],
            s=11 if source != "study" else 20,
            color=color,
            zorder=4,
        )
        ax.text(
            depth[clade] + max_depth * 0.025,
            height[clade],
            label,
            va="center",
            color=color,
            fontsize=10,
            weight="bold" if source == "study" else "normal",
        )
        rows.append(
            {
                "tip": tip,
                "display_label": label,
                "source": source,
                "classification": classification,
                "duplicate_assembly": duplicate,
                "distance_from_subtree_root": depth[clade],
            }
        )
    ax.set(xlim=(-max_depth * 0.055, max_depth * 2.52), ylim=(-1.4, len(tips) - 0.1))
    ax.axis("off")
    scale = 0.1
    ax.plot([0, scale], [-0.9, -0.9], color="#37464A", lw=1.5)
    ax.text(
        scale / 2, -1.12, "0.1 substitutions/site", fontsize=9, ha="center", va="top"
    )
    fig.text(
        0.025,
        0.965,
        "Phylogenomic context of the siphonophore lineages",
        fontsize=12,
        weight="bold",
        va="top",
    )
    fig.text(
        0.025,
        0.916,
        "Metamycoplasmataceae: 21-tip clade from the 483-tip bac120 tree; branch lengths retained",
        fontsize=9,
        color="#48565C",
    )
    handles = [
        Patch(facecolor=color, label=label)
        for color, label in [
            (source_colors["study"], "Siphonophore MAG"),
            (source_colors["external"], "External host-associated genome"),
            (source_colors["GTDB"], "GTDB r220 reference"),
        ]
    ]
    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.065),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    fig.text(
        0.025,
        0.048,
        "Node labels: FastTree local support, not bootstrap percentages. Root support = 0.999.\n"
        "[dup.] The same assembly GCA_929200685.1 occurs as both an external input and a GTDB reference.",
        fontsize=9,
        va="top",
    )
    save_figure(fig, output, "phylogeny_focus")
    write_rows(output / "source_phylogeny_focus.tsv", rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, default=Path("figures/mycoplasmatales"))
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output if args.output.is_absolute() else root / args.output
    output.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42,
                                "font.family": "DejaVu Sans", "font.size": 9})
    inputs = {
        "focal": root / "data/results/mycoplasmatales_gene_content/summary/focal_functions.tsv",
        "genomes": root / "data/results/mycoplasmatales_gene_content/summary/genome_stats.tsv",
        "tree": root / "data/results/mycoplasmatales_phylogeny/mycoplasmatales.bac120.decorated.tree",
        "labels": root / "data/results/mycoplasmatales_phylogeny/tree_labels.tsv",
    }
    before = {str(path.relative_to(root)): checksum(path) for path in inputs.values()}
    genomes = read_rows(inputs["genomes"])
    gene_functions(read_rows(inputs["focal"]), genomes, output)
    focus_tree(inputs["tree"], genomes, read_rows(inputs["labels"]), output)
    if before != {str(path.relative_to(root)): checksum(path) for path in inputs.values()}:
        raise RuntimeError("Source inputs changed during rendering")
    provenance = output / "provenance"
    provenance.mkdir(exist_ok=True)
    snapshot = provenance / Path(__file__).name
    shutil.copyfile(__file__, snapshot)
    artifacts = [path for pattern in ("*.pdf", "*.png", "source_*.tsv")
                 for path in sorted(output.glob(pattern))]
    record = {
        "status": "complete", "run_utc": datetime.now(timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "source_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "source_snapshot": str(snapshot.relative_to(root)), "script_sha256": checksum(snapshot),
        "command": shlex.join([sys.executable, *sys.argv]),
        "environment": {"python": platform.python_version(), "biopython": Bio.__version__,
                        "matplotlib": matplotlib.__version__, "numpy": np.__version__},
        "inputs_sha256": before,
        "outputs_sha256": {str(path.relative_to(root)): checksum(path) for path in artifacts},
        "checks": ["Unique KO/tier records", "Strict assignments are subsets of relaxed assignments",
                   "130 near-complete GTDB reference genomes", "483 unique tree tips; 21 in focal MRCA",
                   "All seven siphonophore MAGs retained", "Input hashes unchanged during rendering"],
        "limits": ["Existing annotations and tree reused, not recomputed",
                   "MAGSP0029 and MAGSP0031 support presence only",
                   "The source tree contains duplicate assembly GCA_929200685.1",
                   "Visual and font inspection recorded separately after rendering"],
    }
    (output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Rendered {len(artifacts)} figures/previews/source tables in {output}")


if __name__ == "__main__":
    main()
