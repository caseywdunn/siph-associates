#!/usr/bin/env python3
"""Plot the full biological scope from accepted tables, without reanalysis.

The microbial display collapses genome targets to explicitly recorded lineages
within each specimen. Eukaryotic cells retain the existing reporting units and
evidence grades. All 205 libraries provide descriptive denominators.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/siph-associates-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle
from matplotlib.transforms import blended_transform_factory

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/overview"
BLUE = "#0072B2"
TEAL = "#009E73"
ORANGE = "#D55E00"
DARK = "#243647"
GRAY = "#EDF1F3"
STUDIES = ["Church2025", "Ahuja2024", "Ahuja2026"]
STUDY_COLORS = {"Church2025": BLUE, "Ahuja2024": "#E69F00", "Ahuja2026": TEAL}
MICROBIAL_GROUPS = [
    "Metamycoplasmataceae",
    "Cognatishimia",
    "Vibrio",
    "CAKMZU01",
    "Pseudoalteromonas",
    "Cutibacterium",
    "Photobacterium",
    "Alteromonas",
    "Remaining bacteria",
    "Archaea",
]
INPUT_NAMES = [
    "figures/cohort/sample_metadata.tsv",
    "figures/cohort/bacterial_detections.tsv",
    "figures/cohort/viral_detections.tsv",
    "figures/cohort/eukaryote_evidence.tsv",
    "figures/collection_context/host_detection_summary.tsv",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(name: str, rows: list[dict[str, object]]) -> Path:
    path = OUT / name
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def save_figure(fig: plt.Figure, name: str) -> list[Path]:
    paths = [OUT / f"{name}.pdf", OUT / f"{name}.png"]
    fig.savefig(paths[0], metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(paths[1], dpi=180)
    plt.close(fig)
    return paths


def host_key(host: str) -> tuple[int, str]:
    return (
        0 if host.startswith("Physalia ") else 1 if host.startswith("Nanomia ") else 2,
        host,
    )


def formatted_host(host: str) -> str:
    """Italicize biological names while retaining provisional qualifiers."""
    if host == "Physonectae":
        return host
    parts = host.split()
    italic_count = 1 if len(parts) == 1 or parts[1].startswith("sp") else 2
    italic = r"\ ".join(parts[:italic_count])
    qualifier = " ".join(parts[italic_count:])
    return rf"$\mathit{{{italic}}}$" + (f" {qualifier}" if qualifier else "")


def panel_title(ax: plt.Axes, letter: str, title: str) -> None:
    ax.set_title(
        f"{letter}  {title}", loc="left", fontsize=11, fontweight="bold", pad=11
    )


def specimen_id(row: dict[str, str]) -> str:
    return row["library_id"].replace(":", "__", 1)


def microbial_group(row: dict[str, str]) -> str:
    parts = row["taxonomy"].split(";")
    if "d__Archaea" in parts:
        return "Archaea"
    if "d__Bacteria" not in parts:
        raise ValueError(f"Unknown bacterial-table domain: {row['taxonomy']}")
    if "f__Metamycoplasmataceae" in parts:
        return "Metamycoplasmataceae"
    genus = next((x[3:] for x in parts if x.startswith("g__")), "")
    return genus if genus in MICROBIAL_GROUPS else "Remaining bacteria"


def draw_overview(
    metadata: list[dict[str, str]], domain_sets: dict[str, set[str]]
) -> list[Path]:
    fig = plt.figure(figsize=(7.8, 8.2))
    grid = fig.add_gridspec(
        3,
        2,
        height_ratios=[1.5, 0.95, 1.1],
        left=0.075,
        right=0.96,
        bottom=0.08,
        top=0.92,
        hspace=0.64,
        wspace=0.4,
    )
    fig.suptitle(
        "Organisms recovered alongside siphonophore genomes",
        fontsize=14,
        fontweight="bold",
        x=0.075,
        ha="left",
        y=0.985,
    )
    photos = fig.add_subplot(grid[0, 0])
    panel_title(photos, "A", "Siphonophore hosts")
    photos.axis("off")
    for x, name, context in [
        (0.0, "Physalia", "Sea surface"),
        (0.54, "Nanomia septata", "Water column"),
    ]:
        photos.add_patch(
            Rectangle(
                (x, 0.25),
                0.46,
                0.7,
                transform=photos.transAxes,
                facecolor="#F7F9FA",
                edgecolor="#8495A1",
                linewidth=1.0,
                linestyle=(0, (4, 3)),
            )
        )
        photos.text(
            x + 0.23,
            0.60,
            "Photograph\nplaceholder",
            ha="center",
            va="center",
            fontsize=10,
            color="#71818C",
            transform=photos.transAxes,
        )
        photos.text(
            x + 0.23,
            0.16,
            name,
            ha="center",
            fontsize=10,
            fontstyle="italic",
            transform=photos.transAxes,
        )
        photos.text(
            x + 0.23,
            0.04,
            context,
            ha="center",
            fontsize=9,
            color=DARK,
            transform=photos.transAxes,
        )

    sampling = fig.add_subplot(grid[0, 1])
    panel_title(sampling, "B", "205 specimens in three studies")
    host_groups = ["Physalia", "Nanomia", "Other siphonophores"]
    sample_rows = []
    left = np.zeros(3)
    for study in STUDIES:
        counts = [
            sum(
                row["study"] == study
                and (
                    row["species_current"].startswith(group + " ")
                    if group != "Other siphonophores"
                    else not row["species_current"].startswith(
                        ("Physalia ", "Nanomia ")
                    )
                )
                for row in metadata
            )
            for group in host_groups
        ]
        sampling.barh(
            range(3), counts, left=left, color=STUDY_COLORS[study], height=0.55
        )
        for group, count in zip(host_groups, counts):
            sample_rows.append(
                {"host_group": group, "study": study, "libraries": count}
            )
        left += counts
    for y, count in enumerate(left):
        sampling.text(count + 3, y, str(int(count)), va="center", fontsize=10)
    sampling.set_yticks(
        range(3), ["Physalia", "Nanomia", "Other\nsiphonophores"], fontsize=10
    )
    sampling.invert_yaxis()
    sampling.set_xlim(0, 170)
    sampling.set_xticks([0, 50, 100, 150])
    sampling.set_xlabel("Specimens", fontsize=10)
    sampling.spines[["left", "right", "top"]].set_visible(False)
    sampling.tick_params(axis="y", length=0)
    sampling.legend(
        handles=[Patch(color=STUDY_COLORS[s], label=s) for s in STUDIES],
        frameon=False,
        loc="upper right",
        fontsize=9,
        bbox_to_anchor=(1.05, -0.23),
        ncol=1,
    )

    workflow = fig.add_subplot(grid[1, :])
    panel_title(workflow, "C", "Host sequencing provides a second biological record")
    workflow.axis("off")
    labels = [
        "Collected\nsiphonophore\ntissue",
        "Short-read\nhost-genome\nlibraries",
        "Identify and\nvalidate non-host\nsequences",
        "Interpret organisms\nusing natural\nhistory",
    ]
    for index, label in enumerate(labels):
        x = 0.01 + index * 0.257
        workflow.add_patch(
            FancyBboxPatch(
                (x, 0.10),
                0.21,
                0.72,
                boxstyle="round,pad=0.012",
                facecolor=GRAY,
                edgecolor="none",
                transform=workflow.transAxes,
            )
        )
        workflow.text(
            x + 0.105,
            0.46,
            label,
            ha="center",
            va="center",
            fontsize=10,
            transform=workflow.transAxes,
            linespacing=1.4,
        )
        if index < 3:
            workflow.add_patch(
                FancyArrowPatch(
                    (x + 0.222, 0.46),
                    (x + 0.25, 0.46),
                    arrowstyle="-|>",
                    mutation_scale=12,
                    color="#5C6D78",
                    transform=workflow.transAxes,
                )
            )

    detection = fig.add_subplot(grid[2, 0])
    panel_title(detection, "D", "Recovered across the cohort")
    domain_colors = [BLUE, "#9B6BAC", TEAL, ORANGE]
    for index, (domain, samples) in enumerate(domain_sets.items()):
        detection.barh(index, 205, color=GRAY, height=0.55)
        detection.barh(index, len(samples), color=domain_colors[index], height=0.55)
        detection.text(211, index, f"{len(samples)}/205", va="center", fontsize=10)
    detection.set_yticks(range(4), list(domain_sets), fontsize=10)
    detection.invert_yaxis()
    detection.set_xlim(0, 270)
    detection.set_xticks([0, 100, 205])
    detection.set_xlabel("Specimens with supported detections", fontsize=9)
    detection.tick_params(axis="y", length=0)
    detection.spines[["left", "right", "top"]].set_visible(False)
    position = detection.get_position()
    detection.set_position(
        [position.x0 + 0.055, position.y0, position.width - 0.055, position.height]
    )

    questions = fig.add_subplot(grid[2, 1])
    panel_title(questions, "E", "Biological questions")
    questions.axis("off")
    for y, title, description in [
        (
            0.88,
            "Microbial associates",
            "Which bacteria, archaea and viruses\naccompany siphonophores?",
        ),
        (
            0.48,
            "Prey and parasites",
            "Which eukaryotes suggest feeding\nor parasite associations?",
        ),
        (
            0.08,
            "Host and collection context",
            "How do detections differ among\nhosts and sampling locations?",
        ),
    ]:
        questions.text(0, y, title, fontsize=10, fontweight="bold", color=DARK)
        questions.text(0, y - 0.08, description, va="top", fontsize=9, linespacing=1.3)
    return save_figure(fig, "study_overview") + [
        write_tsv("study_sampling.tsv", sample_rows)
    ]


def draw_microbes(
    metadata: list[dict[str, str]], bacteria: list[dict[str, str]]
) -> list[Path]:
    meta = {specimen_id(row): row for row in metadata}
    host_counts = Counter(row["species_current"] for row in metadata)
    hosts = sorted(host_counts, key=host_key)
    present: dict[tuple[str, str], set[str]] = defaultdict(set)
    any_bacteria: dict[str, set[str]] = defaultdict(set)
    classified_rows = []
    for row in bacteria:
        sample = row["sample_id"]
        host = meta[sample]["species_current"]
        group = microbial_group(row)
        present[(host, group)].add(sample)
        if group != "Archaea":
            any_bacteria[host].add(sample)
        classified_rows.append(
            {
                "sample_id": sample,
                "host": host,
                "target_id": row["target_id"],
                "display_group": group,
                "taxonomy": row["taxonomy"],
            }
        )
    counts = np.array(
        [[len(present[(host, group)]) for group in MICROBIAL_GROUPS] for host in hosts]
    )
    denoms = np.array([host_counts[host] for host in hosts])
    fig = plt.figure(figsize=(7.8, 8.2))
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[10, 1.3],
        left=0.285,
        right=0.955,
        bottom=0.165,
        top=0.76,
        wspace=0.12,
    )
    ax = fig.add_subplot(grid[0, 0])
    totals = fig.add_subplot(grid[0, 1], sharey=ax)
    im = ax.imshow(
        counts / denoms[:, None], aspect="auto", vmin=0, vmax=1, cmap="Blues"
    )
    for i in range(len(hosts)):
        for j in range(len(MICROBIAL_GROUPS)):
            if counts[i, j]:
                ax.text(
                    j,
                    i,
                    str(counts[i, j]),
                    ha="center",
                    va="center",
                    fontsize=9,
                    color="white" if counts[i, j] / denoms[i] > 0.55 else DARK,
                )
    ax.set_yticks(
        range(len(hosts)),
        [f"{formatted_host(h)} ({host_counts[h]})" for h in hosts],
        fontsize=9,
    )
    ax.set_xticks(
        range(len(MICROBIAL_GROUPS)),
        MICROBIAL_GROUPS,
        rotation=57,
        ha="left",
        rotation_mode="anchor",
        fontsize=9,
    )
    ax.xaxis.tick_top()
    for group, label in zip(MICROBIAL_GROUPS, ax.get_xticklabels()):
        if group not in {"Metamycoplasmataceae", "Remaining bacteria", "Archaea"}:
            label.set_fontstyle("italic")
    ax.tick_params(axis="both", length=0, pad=3)
    ax.set_xticks(np.arange(-0.5, len(MICROBIAL_GROUPS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(hosts), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False, top=False)
    ax.spines[:].set_visible(False)
    ax.axvline(8.5, color=DARK, linewidth=1.1)
    totals.barh(
        range(len(hosts)),
        [len(any_bacteria[h]) / host_counts[h] for h in hosts],
        height=0.65,
        color="#687F91",
    )
    totals.set_xlim(0, 1.1)
    totals.set_xticks([0, 1], ["0", "1"])
    totals.tick_params(axis="y", left=False, labelleft=False)
    totals.spines[["top", "right", "left"]].set_visible(False)
    totals.text(
        0.5,
        1.02,
        "Any\nbacteria",
        transform=totals.transAxes,
        ha="center",
        va="bottom",
        fontsize=9,
    )
    for y, host in enumerate(hosts):
        if y and host_key(host)[0] != host_key(hosts[y - 1])[0]:
            ax.axhline(y - 0.5, color=DARK, linewidth=1)
            totals.axhline(y - 0.5, color=DARK, linewidth=1)
    fig.suptitle(
        "Bacteria and archaea across sampled siphonophores",
        fontsize=13,
        fontweight="bold",
        x=0.055,
        ha="left",
        y=0.975,
    )
    fig.text(
        0.055,
        0.929,
        "All 205 specimens retained • host sample sizes in parentheses",
        fontsize=9,
    )
    color_ax = fig.add_axes([0.31, 0.12, 0.37, 0.017])
    cbar = fig.colorbar(
        im, cax=color_ax, orientation="horizontal", ticks=[0, 0.25, 0.5, 0.75, 1]
    )
    cbar.set_label("Fraction of specimens with a supported detection", fontsize=9)
    cbar.ax.tick_params(labelsize=9)
    fig.text(
        0.055,
        0.025,
        "Numbers inside cells: positive specimens. Empty cells: no retained detection.\n"
        "Selected bacterial lineages are shown separately; remaining bacteria are pooled.",
        fontsize=9,
        linespacing=1.45,
    )
    data = [
        {
            "host": host,
            "group": group,
            "positive_libraries": int(counts[i, j]),
            "total_libraries": host_counts[host],
            "fraction_positive": counts[i, j] / host_counts[host],
        }
        for i, host in enumerate(hosts)
        for j, group in enumerate(MICROBIAL_GROUPS)
    ]
    return save_figure(fig, "microbial_host_distribution") + [
        write_tsv("microbial_host_distribution.tsv", data),
        write_tsv("microbial_display_assignments.tsv", classified_rows),
    ]


def draw_eukaryotes(
    metadata: list[dict[str, str]], retained: list[dict[str, str]]
) -> list[Path]:
    meta = {specimen_id(row): row for row in metadata}
    units = sorted({r["reporting_unit"] for r in retained})
    samples = sorted(
        {r["sample_id"] for r in retained},
        key=lambda s: (host_key(meta[s]["species_current"]), meta[s]["study"], s),
    )
    unit_index = {unit: i for i, unit in enumerate(units)}
    sample_index = {sample: i for i, sample in enumerate(samples)}
    codes = {"validated": 1, "high_confidence": 2}
    matrix = np.zeros((len(samples), len(units)), dtype=int)
    display_rows = []
    for row in retained:
        i, j = sample_index[row["sample_id"]], unit_index[row["reporting_unit"]]
        if matrix[i, j]:
            raise ValueError("Duplicate retained eukaryotic unit within a library")
        matrix[i, j] = codes[row["grade"]]
        display_rows.append(
            {
                "sample_id": row["sample_id"],
                "host": meta[row["sample_id"]]["species_current"],
                "display_row": i + 1,
                "reporting_unit": row["reporting_unit"],
                "grade": row["grade"],
            }
        )
    if (int(np.count_nonzero(matrix)), len(samples), int((matrix == 2).sum())) != (
        61,
        39,
        42,
    ):
        raise ValueError(
            "Review plot scope: expected 61 detections/39 libraries/42 read+assembly cells"
        )
    labels = [unit.split(";")[-1] for unit in units]
    if len(set(labels)) != len(labels):
        raise ValueError("Ambiguous shortened eukaryotic labels")
    fig = plt.figure(figsize=(10.4, 10.4))
    grid = fig.add_gridspec(
        1,
        2,
        width_ratios=[16, 2.2],
        left=0.36,
        right=0.955,
        bottom=0.08,
        top=0.705,
        wspace=0.10,
    )
    ax = fig.add_subplot(grid[0, 0])
    total = fig.add_subplot(grid[0, 1], sharey=ax)
    ax.imshow(
        matrix,
        aspect="auto",
        cmap=ListedColormap(["#F2F3F4", "#56B4E9", BLUE]),
        norm=BoundaryNorm([-0.5, 0.5, 1.5, 2.5], 3),
    )
    ax.set_yticks(
        range(len(samples)),
        [meta[s]["library_id"].split(":", 1)[1] for s in samples],
        fontsize=12,
    )
    ax.set_xticks(
        range(len(units)),
        labels,
        rotation=57,
        ha="left",
        rotation_mode="anchor",
        fontsize=12,
    )
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=3)
    ax.tick_params(axis="y", length=0, pad=14)
    ax.set_xticks(np.arange(-0.5, len(units), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(samples), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.6)
    ax.tick_params(which="minor", bottom=False, left=False, top=False)
    ax.spines[:].set_visible(False)
    host_groups: dict[str, list[int]] = defaultdict(list)
    for i, sample in enumerate(samples):
        host_groups[meta[sample]["species_current"]].append(i)
        ax.plot(
            -0.78,
            i,
            marker="s",
            markersize=5.8,
            color=STUDY_COLORS[meta[sample]["study"]],
            clip_on=False,
        )
    for host, indices in host_groups.items():
        ax.text(
            0.035,
            (indices[0] + indices[-1]) / 2,
            formatted_host(host),
            ha="left",
            va="center",
            fontsize=12,
            transform=blended_transform_factory(fig.transFigure, ax.transData),
        )
        if indices[0]:
            ax.axhline(indices[0] - 0.5, color=DARK, linewidth=1.0)
            total.axhline(indices[0] - 0.5, color=DARK, linewidth=1.0)
    # Markers at the left of the matrix must not expand the data limits.
    ax.set_xlim(-0.5, len(units) - 0.5)
    ax.set_ylim(len(samples) - 0.5, -0.5)
    total.barh(
        range(len(samples)), (matrix == 1).sum(axis=1), color="#56B4E9", height=0.70
    )
    total.barh(
        range(len(samples)),
        (matrix == 2).sum(axis=1),
        left=(matrix == 1).sum(axis=1),
        color=BLUE,
        height=0.70,
    )
    total.set_xlim(0, 5.4)
    total.set_xticks([0, 2, 4])
    total.set_xlabel("Groups", fontsize=12)
    total.tick_params(axis="x", labelsize=12)
    total.tick_params(axis="y", left=False, labelleft=False)
    total.spines[["top", "right", "left"]].set_visible(False)
    total.text(
        0.5,
        1.02,
        "Groups per\nlibrary",
        transform=total.transAxes,
        ha="center",
        va="bottom",
        fontsize=12,
    )
    fig.suptitle(
        "Eukaryotic sequences across siphonophore hosts",
        x=0.035,
        ha="left",
        y=0.975,
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.035,
        0.933,
        "61 eukaryote group detections in 39 of 205 libraries",
        fontsize=12,
    )
    fig.legend(
        handles=[
            Patch(color="#56B4E9", label="Read or assembly support"),
            Patch(color=BLUE, label="Read and assembly support"),
        ],
        loc="upper left",
        bbox_to_anchor=(0.028, 0.905),
        ncol=2,
        frameon=False,
        fontsize=12,
    )
    fig.legend(
        handles=[Patch(color=STUDY_COLORS[s], label=s) for s in STUDIES],
        loc="upper left",
        bbox_to_anchor=(0.028, 0.867),
        ncol=3,
        frameon=False,
        fontsize=12,
    )
    fig.text(
        0.035,
        0.025,
        "Rows show the 39 libraries with retained eukaryotic evidence, grouped by host.\n"
        "Grey cells indicate no retained detection; broad groups may contain more than one organism.",
        fontsize=12,
        linespacing=1.5,
    )
    return save_figure(fig, "eukaryote_evidence") + [
        write_tsv("eukaryote_display.tsv", display_rows)
    ]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    tables = {Path(name).name: read_tsv(ROOT / name) for name in INPUT_NAMES}
    metadata = tables["sample_metadata.tsv"]
    meta = {specimen_id(row): row for row in metadata}
    if len(metadata) != 205 or len(meta) != 205:
        raise ValueError("Expected one metadata row for each of 205 libraries")
    bacteria = tables["bacterial_detections.tsv"]
    viruses = tables["viral_detections.tsv"]
    eukaryotes = [
        row
        for row in tables["eukaryote_evidence.tsv"]
        if row["count_as_detection"] == "true"
    ]
    for rows in [bacteria, viruses, eukaryotes]:
        for row in rows:
            if row["sample_id"] not in meta or row["grade"] not in {
                "validated",
                "high_confidence",
            }:
                raise ValueError(f"Unexpected sample or evidence grade: {row}")
    domain_sets = {
        "Bacteria": {
            r["sample_id"]
            for r in bacteria
            if "d__Bacteria" in r["taxonomy"].split(";")
        },
        "Archaea": {
            r["sample_id"] for r in bacteria if "d__Archaea" in r["taxonomy"].split(";")
        },
        "Viruses": {r["sample_id"] for r in viruses},
        "Eukaryotes": {r["sample_id"] for r in eukaryotes},
    }
    # Check the new display's independent collapse against accepted host totals.
    for field, name in [
        ("bacterial", "Bacteria"),
        ("archaeal", "Archaea"),
        ("viral", "Viruses"),
        ("eukaryote", "Eukaryotes"),
    ]:
        count = sum(
            int(r[f"{field}_positive_libraries"])
            for r in tables["host_detection_summary.tsv"]
        )
        if count != len(domain_sets[name]):
            raise ValueError(
                f"Domain summary disagrees with accepted host table: {name}"
            )
    outputs = draw_overview(metadata, domain_sets)
    outputs += draw_microbes(metadata, bacteria)
    outputs += draw_eukaryotes(metadata, eukaryotes)
    outputs.append(
        write_tsv(
            "domain_summary.tsv",
            [
                {
                    "domain": domain,
                    "positive_libraries": len(samples),
                    "total_libraries": 205,
                }
                for domain, samples in domain_sets.items()
            ],
        )
    )
    script = Path(__file__).resolve()
    snapshot = OUT / "provenance/plot_associate_overview.py"
    snapshot.parent.mkdir(exist_ok=True)
    shutil.copyfile(script, snapshot)
    manifest = {
        "status": "complete",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "command": shlex.join([sys.executable, *sys.argv]),
        "working_directory": str(ROOT),
        "source_snapshot": str(snapshot.relative_to(ROOT)),
        "script_sha256": sha256(script),
        "inputs_sha256": {name: sha256(ROOT / name) for name in INPUT_NAMES},
        "outputs_sha256": {
            str(path.relative_to(ROOT)): sha256(path) for path in outputs
        },
        "environment": {
            "python": platform.python_version(),
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
        },
        "scope": "All 205 libraries retained descriptively; no flowcell-based inference or sequence reanalysis.",
        "checks": {
            "metadata_libraries": 205,
            "eukaryotic_detections": 61,
            "eukaryote_positive_libraries": 39,
            "eukaryote_read_and_assembly_cells": 42,
            "domain_counts_match_accepted_host_summary": True,
        },
        "microbial_display": "Mutually exclusive target assignments: Metamycoplasmataceae, seven frequent named bacterial groups, remaining bacteria, and archaea; presence deduplicated within each specimen.",
        "photographs": "Explicit empty placeholders; no photographs or specimen images supplied.",
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["checks"], indent=2))


if __name__ == "__main__":
    main()
