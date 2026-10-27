#!/usr/bin/env python3
"""Render publication figures from retained Mycoplasmatales results.

Run from the repository root on a compute node. This does not infer a new tree,
or change genome classifications. Source tables and input checksums
are retained beside the figures.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
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
from matplotlib.patches import Patch

MAG_ORDER = [
    "MAGSP0005", "MAGSP0007", "MAGSP0010", "MAGSP0029",
    "MAGSP0031", "MAGSP0011", "MAGSP0012",
]
IQTREE_MODEL = "LG+F+I+R10"
WEAK_RELATIONSHIP_TIPS = frozenset({
    "MAGSP0007", "MAGSP0010", "MAGSP0029", "MAGSP0031", "MAGSP0011",
    "GB_GCA_013215685.1", "GB_GCA_929200685.1", "EXT_GCA_060113895.1",
    "EXT_GCA_929200505.1", "EXT_GCA_929200285.1", "EXT_GCA_929200685.1",
})


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


def node_support(clade) -> tuple[float, float] | None:
    """Read IQ-TREE's SH-aLRT/UFBoot percentages without conflating the measures."""
    if clade.is_terminal():
        return None
    if clade.name is None and clade.confidence is None:
        return None
    fields = (clade.name or str(clade.confidence)).split("/")
    if len(fields) != 2:
        raise ValueError(f"Expected SH-aLRT/UFBoot support: {clade.name}")
    values = tuple(float(value) for value in fields)
    if any(not math.isfinite(value) or not 0 <= value <= 100 for value in values):
        raise ValueError(f"Invalid support percentages: {clade.name}")
    return values


def split_key(tips: frozenset[str], universe: frozenset[str]) -> tuple[str, ...]:
    """Identify an edge independently of the root and the chosen descendant side."""
    sides = (tuple(sorted(tips)), tuple(sorted(universe - tips)))
    return min(sides, key=lambda side: (len(side), side))


def root_on_outgroup(tree, outgroup: frozenset[str]) -> None:
    """Root at the outgroup edge, preserving support and lengths on the same splits."""
    tips = [tip.name for tip in tree.get_terminals()]
    universe = frozenset(tips)
    if len(tips) != len(universe) or not outgroup or not outgroup < universe:
        raise ValueError("Tree tips must be unique and contain the complete outgroup")
    edge_support, edge_lengths = {}, {}
    target = None
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        descendants = frozenset(tip.name for tip in clade.get_terminals())
        key = split_key(descendants, universe)
        edge_lengths[key] = edge_lengths.get(key, 0.0) + (clade.branch_length or 0.0)
        support = node_support(clade)
        if support is not None:
            if key in edge_support and edge_support[key] != support:
                raise ValueError(f"Conflicting support on an equivalent edge: {key}")
            edge_support[key] = support
        if descendants in (outgroup, universe - outgroup):
            target = clade
    if target is None:
        raise ValueError("The complete outgroup does not define a split in this tree")

    # IQ-TREE can display the unrooted tree at one outgroup tip. Rerooting changes
    # parent/child directions, so retain support by edge identity, not node object.
    tree.root_with_outgroup(target, outgroup_branch_length=0.0)
    tree.root.name = None
    tree.root.confidence = None
    after_lengths = {}
    for clade in tree.find_clades():
        if clade is tree.root:
            continue
        descendants = frozenset(tip.name for tip in clade.get_terminals())
        key = split_key(descendants, universe)
        after_lengths[key] = after_lengths.get(key, 0.0) + (clade.branch_length or 0.0)
        if not clade.is_terminal():
            support = edge_support.get(key)
            clade.name = "/".join(f"{value:g}" for value in support) if support else None
            clade.confidence = None
    if edge_lengths.keys() != after_lengths.keys() or any(
        not math.isclose(length, after_lengths[key], rel_tol=1e-10, abs_tol=1e-12)
        for key, length in edge_lengths.items()
    ):
        raise ValueError("Rerooting changed an unrooted edge or its length")
    root_sides = {frozenset(tip.name for tip in child.get_terminals())
                  for child in tree.root.clades}
    if root_sides != {outgroup, universe - outgroup}:
        raise ValueError("Rerooting did not separate the entire outgroup")


def source_host(
    tip: str, external: dict[str, dict], samples: dict[str, dict],
    reviewed_hosts: dict[str, dict],
) -> dict[str, str]:
    """Join specimen or reviewed accession records, including environmental sources."""
    if tip.startswith("MAGSP"):
        source_record = external[tip]["host_or_source"]
        if not source_record.startswith("source library "):
            raise ValueError(f"Missing source-library identifier for {tip}")
        library = source_record.removeprefix("source library ")
        if library not in samples or not samples[library]["species_current"]:
            raise ValueError(f"Missing source-specimen host for {tip}: {library}")
        return {
            "source_host": samples[library]["species_current"],
            "host_clade": "Siphonophora",
            "source_type": "host",
            "source_host_basis": "source specimen metadata",
            "source_host_record": library,
            "source_host_url": "figures/cohort/sample_metadata.tsv",
            "host_clade_basis": "siphonophore source specimen identification",
        }
    accession = tip.removeprefix("EXT_").removeprefix("RS_").removeprefix("GB_")
    if accession not in reviewed_hosts:
        raise ValueError(f"Missing reviewed host/source record for {accession}")
    return {
        key: reviewed_hosts[accession][key]
        for key in ("source_host", "host_clade", "source_type", "source_host_basis",
                    "source_host_record", "source_host_url", "host_clade_basis")
    }


def focus_tree(
    tree_path: Path, outgroup_path: Path, genomes: list[dict], labels: list[dict],
    sample_metadata: list[dict], host_records: list[dict], output: Path
) -> None:
    """Display the existing MRCA subtree without removing duplicate source tips."""
    tree = Phylo.read(tree_path, "newick")
    outgroup = frozenset(outgroup_path.read_text().strip().split(","))
    if len(outgroup) != 94:
        raise ValueError("Expected the preserved 94-genome Acholeplasmataceae outgroup")
    root_on_outgroup(tree, outgroup)
    terminals = tree.get_terminals()
    names = [tip.name for tip in terminals]
    if len(names) != len(set(names)) or not set(MAG_ORDER).issubset(names):
        raise ValueError("Tree tips are duplicated or missing selected study MAGs")
    root = tree.common_ancestor(MAG_ORDER)
    root.ladderize()
    tips = root.get_terminals()
    if len(terminals) != 483 or len(tips) != 21:
        raise ValueError(
            "Tree changed: inspect the selected 21-tip clade before plotting"
        )
    taxonomy = {row["genome"]: row for row in [*genomes, *labels]}
    external = {row["genome"]: row for row in labels}
    samples = {row["analysis_sample_id"]: row for row in sample_metadata}
    if len(samples) != len(sample_metadata):
        raise ValueError("Source specimen metadata contain duplicate analysis sample IDs")
    reviewed_hosts = {row["accession"]: row for row in host_records}
    if len(reviewed_hosts) != len(host_records):
        raise ValueError("Reviewed source metadata contain duplicate assembly accessions")
    for accession, row in reviewed_hosts.items():
        if not all(row[key] for key in ("source_host", "host_clade", "source_type",
                                        "source_host_basis", "source_host_record",
                                        "source_host_url", "host_clade_basis")):
            raise ValueError(f"Incomplete reviewed source record: {accession}")
        if row["source_type"] not in {"host", "environmental", "unresolved"}:
            raise ValueError(f"Unrecognized source type for {accession}")
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
    fig, (ax, label_ax, host_ax, clade_ax) = plt.subplots(
        ncols=4, sharey=True, figsize=(11.5, 8.2),
        gridspec_kw={"width_ratios": [2.7, 3.6, 3.2, 1.7], "wspace": 0.015},
    )
    fig.subplots_adjust(left=0.045, right=0.985, bottom=0.14, top=0.92)
    rows, support_rows = [], []
    source_colors = {"study": "#197F8C", "external": "#B75B24", "GTDB": "#48565C"}
    mag_lineages = {
        "MAGSP0005": "clade A",
        "MAGSP0007": "clade B",
        "MAGSP0010": "clade B",
        "MAGSP0029": "clade B",
        "MAGSP0031": "sister to clade B",
        "MAGSP0011": "DT-68",
        "MAGSP0012": "Mycoplasma_K",
    }
    for clade in root.find_clades(order="preorder"):
        if clade.clades:
            ax.plot(
                [depth[clade]] * 2,
                [height[clade.clades[0]], height[clade.clades[-1]]],
                color="black",
                linewidth=0.8,
            )
            for child in clade.clades:
                ax.plot(
                    [depth[clade], depth[child]],
                    [height[child]] * 2,
                    color="black",
                    linewidth=0.8,
                )
            support = node_support(clade)
            if support is not None:
                support_rows.append({
                    "node_id": f"node{len(support_rows) + 1:02d}",
                    "descendant_tips": ";".join(sorted(tip.name for tip in clade.get_terminals())),
                    "tip_count": len(clade.get_terminals()),
                    "sh_alrt_percent": support[0],
                    "ultrafast_bootstrap_percent": support[1],
                    "branch_length": clade.branch_length,
                    "distance_from_subtree_root": depth[clade],
                    "display_y": height[clade],
                })
                if support == (100.0, 100.0):
                    ax.scatter(depth[clade], height[clade], s=12, color="black", zorder=5)
                else:
                    ax.text(
                        depth[clade] - 0.008 * max_depth,
                        height[clade] + 0.16,
                        f"{support[0]:g}/{support[1]:g}",
                        fontsize=10,
                        ha="right",
                        va="bottom",
                        color="#37464A",
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
            label = tip.removeprefix("EXT_")
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
        host = source_host(tip, external, samples, reviewed_hosts)
        color = source_colors[source]
        # Shared y coordinates keep each host entry on its corresponding tip row.
        ax.plot(
            [depth[clade], max_depth * 1.05], [height[clade]] * 2,
            color="#D8DEE0", linewidth=0.45, linestyle=(0, (1.5, 2.5)), zorder=0,
        )
        ax.scatter(
            depth[clade],
            height[clade],
            s=11 if source != "study" else 20,
            color=color,
            zorder=4,
        )
        label_ax.text(
            0.02,
            height[clade],
            label,
            va="center",
            color=color,
            fontsize=11,
            weight="bold" if source == "study" else "normal",
        )
        host_ax.text(
            0.045, height[clade], host["source_host"], va="center",
            color="#37464A", fontsize=11,
            style="italic" if host["source_type"] == "host" else "normal",
        )
        clade_ax.text(
            0.045, height[clade], host["host_clade"], va="center",
            color="#37464A", fontsize=11,
        )
        rows.append(
            {
                "tip": tip,
                "display_label": label,
                "source": source,
                "classification": classification,
                "duplicate_assembly": duplicate,
                "inference_method": "IQ-TREE 3.0.1",
                "substitution_model": IQTREE_MODEL,
                "distance_from_subtree_root": depth[clade],
                "display_y": height[clade],
                **host,
            }
        )
    ax.set(xlim=(-max_depth * 0.24, max_depth * 1.05), ylim=(-1.4, len(tips) + 0.8))
    for table_ax, heading in ((label_ax, "Genome"), (host_ax, "Source host / habitat"),
                              (clade_ax, "Host clade")):
        table_ax.set_xlim(0, 1)
        table_ax.text(0.02 if table_ax is label_ax else 0.045, len(tips) + 0.05,
                      heading, fontsize=11, weight="bold", color="#37464A",
                      va="center")
        table_ax.axhline(len(tips) - 0.55, color="#65767B", linewidth=0.7)
        table_ax.axhline(-0.5, color="#65767B", linewidth=0.7)
        for row in rows:
            if row["display_y"] % 2 == 0:
                table_ax.axhspan(row["display_y"] - 0.5, row["display_y"] + 0.5,
                                color="#F5F7F7", zorder=0)
    for panel in (ax, label_ax, host_ax, clade_ax):
        panel.axis("off")
    scale = 0.1
    ax.plot([0, scale], [-0.9, -0.9], color="#37464A", lw=1.5)
    ax.text(
        scale / 2, -1.12, "0.1 substitutions/site", fontsize=9, ha="center", va="top"
    )
    handles = [
        Patch(facecolor=color, label=label)
        for color, label in [
            (source_colors["study"], "Bacterial MAG (this study)"),
            (source_colors["external"], "External host-associated genome"),
            (source_colors["GTDB"], "GTDB r220 reference"),
        ]
    ]
    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.025, 0.035),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    save_figure(fig, output, "phylogeny_focus")
    write_rows(output / "source_phylogeny_focus.tsv", rows)
    write_rows(output / "source_phylogeny_support.tsv", support_rows)
    write_rows(output / "source_phylogeny_hosts.tsv", [
        {key: row[key] for key in ("tip", "display_y", "source", "source_host",
                                  "host_clade", "source_type", "source_host_basis",
                                  "source_host_record", "source_host_url", "host_clade_basis")}
        for row in rows
    ])
    weak_clade = tree.common_ancestor(sorted(WEAK_RELATIONSHIP_TIPS))
    if frozenset(tip.name for tip in weak_clade.get_terminals()) != WEAK_RELATIONSHIP_TIPS:
        raise ValueError("The clade-B/MAGSP0031/DT-68 relationship changed")
    focal_support, weak_support = node_support(root), node_support(weak_clade)
    if focal_support is None or weak_support is None:
        raise ValueError("Support is missing for a relationship reported in the manuscript")
    values = {
        "phylogenyFocalSHaLRT": focal_support[0],
        "phylogenyFocalBootstrap": focal_support[1],
        "phylogenyWeakSHaLRT": weak_support[0],
        "phylogenyWeakBootstrap": weak_support[1],
    }
    (output / "values.tex").write_text(
        "% Generated from the accepted IQ-TREE tree by plot_mycoplasmatales_figures.py.\n"
        + "".join(f"\\newcommand{{\\{key}}}{{{value:g}}}\n" for key, value in values.items())
    )


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
        "genomes": root / "config/mycoplasmatales_reference_taxonomy.tsv",
        "tree": root / "data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/bac120.treefile",
        "outgroup": root / "data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/input/outgroup.txt",
        "inference_report": root / "data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/bac120.iqtree",
        "inference_complete": root / "data/results/mycoplasmatales_phylogeny/iqtree3/2026-10-25/inference_complete.json",
        "labels": root / "data/results/mycoplasmatales_phylogeny/tree_labels.tsv",
        "sample_metadata": root / "figures/cohort/sample_metadata.tsv",
        "host_records": root / "config/mycoplasmatales_source_hosts.tsv",
    }
    before = {str(path.relative_to(root)): checksum(path) for path in inputs.values()}
    completion = json.loads(inputs["inference_complete"].read_text())
    if completion["tips"] != 483 or completion["independent_searches"] != 5:
        raise ValueError("The accepted inference must contain 483 tips and five searches")
    if f"Model of substitution: {IQTREE_MODEL}" not in inputs["inference_report"].read_text():
        raise ValueError("The inferred substitution model differs from the reviewed model")
    genomes = read_rows(inputs["genomes"])
    focus_tree(inputs["tree"], inputs["outgroup"], genomes, read_rows(inputs["labels"]),
               read_rows(inputs["sample_metadata"]), read_rows(inputs["host_records"]), output)
    if before != {str(path.relative_to(root)): checksum(path) for path in inputs.values()}:
        raise RuntimeError("Source inputs changed during rendering")
    provenance = output / "provenance"
    provenance.mkdir(exist_ok=True)
    snapshot = provenance / Path(__file__).name
    shutil.copyfile(__file__, snapshot)
    artifacts = [output / name for name in ("phylogeny_focus.pdf",
                 "phylogeny_focus.png", "source_phylogeny_focus.tsv", "source_phylogeny_hosts.tsv",
                 "source_phylogeny_support.tsv", "values.tex")]
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
        "phylogeny": {"software": "IQ-TREE 3.0.1", "model": IQTREE_MODEL,
                      "support": ["SH-aLRT percent", "ultrafast bootstrap percent"],
                      "outgroup": "94 Acholeplasmataceae genomes", "independent_searches": 5,
                      "support_replicates_per_search": 1000,
                      "retained_result": "Highest-likelihood search; supports are not pooled"},
        "outputs_sha256": {str(path.relative_to(root)): checksum(path) for path in artifacts},
        "checks": ["483 unique tree tips; 21 in focal MRCA",
                   "94-tip outgroup split verified; edge lengths and support preserved on rerooting",
                   "All seven siphonophore MAGs retained", "Input hashes unchanged during rendering",
                   "Source-host and clade rows share exact tree-tip y coordinates",
                   "Source hosts joined by source-library ID or identical assembly accession",
                   "Every reference accession has a reviewed host or environmental source record"],
        "limits": ["Completed IQ-TREE inference and GTDB r220 classifications reused, not recomputed",
                   "The source tree contains duplicate assembly GCA_929200685.1",
                   "Environmental samples do not establish an organismal host; host clade shown as a dash",
                   "Visual and font inspection recorded separately after rendering"],
    }
    (output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n")
    print(f"Rendered {len(artifacts)} figures/previews/source tables in {output}")


if __name__ == "__main__":
    main()
