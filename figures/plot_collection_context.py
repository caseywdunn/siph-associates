#!/usr/bin/env python3
"""Plot retained collection-context tables without rerunning sequence analyses.

All biological detections and grouping decisions are read from the companion
analysis tables. Plotting only arranges those values, labels, and denominators.
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

os.environ.setdefault("MPLCONFIGDIR", "/tmp/siph-associates-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
BLUE = "#0072B2"
DARK = "#243647"
ORANGE = "#D55E00"
PALE = "#F1F5F7"
DOMAIN_FIELDS = (
    ("bacterial", "Bacteria"),
    ("archaeal", "Archaea"),
    ("viral", "Viruses"),
    ("eukaryote", "Eukaryotes"),
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a retained table; malformed or empty inputs should fail visibly."""
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No rows in required plot input: {path}")
    return rows


def number(value: str) -> float | None:
    return float(value) if value not in ("", "NA", "NaN", "None") else None


def is_true(value: str) -> bool:
    return value.lower() in ("true", "1", "yes")


def abbreviated_sample(row: dict[str, str]) -> str:
    return row["library_id"].split(":", 1)[-1]


def depth_sort(row: dict[str, str]) -> tuple[float, str]:
    value = number(row["depth_m"])
    if value is None:
        value = number(row["depth_min_m"])
    return (float("inf") if value is None else value, row["library_id"])


def depth_axis(ax: plt.Axes, maximum: float) -> None:
    """Keep surface and shallow collection intervals legible on the same axis."""
    ax.set_xscale("symlog", linthresh=20, linscale=0.7)
    ticks = [x for x in (0, 20, 100, 1000) if x <= maximum]
    ax.set_xticks(ticks, [str(x) for x in ticks])
    ax.set_xlim(-2, maximum)
    ax.set_xlabel("Collection depth (m; linear to 20 m, then logarithmic)", fontsize=9)
    ax.grid(axis="x", color="#DDE3E7", linewidth=0.7, zorder=0)
    ax.tick_params(axis="y", length=0)
    ax.spines[["top", "right", "left"]].set_visible(False)


def draw_depth(
    ax: plt.Axes,
    unknown_ax: plt.Axes,
    row: dict[str, str],
    y: float,
    color: str = BLUE,
    filled: bool = True,
    marker: str = "o",
    size: float = 4.2,
) -> None:
    point = number(row["depth_m"])
    low = number(row["depth_min_m"])
    high = number(row["depth_max_m"])
    if low is not None and high is not None:
        ax.plot([low, high], [y, y], color=color, linewidth=1.5, zorder=3)
        ax.plot(
            [low, high], [y, y], linestyle="none", marker="|", color=color, markersize=5
        )
    if point is not None:
        ax.plot(
            point,
            y,
            marker=marker,
            markersize=size,
            markerfacecolor=color if filled else "white",
            markeredgecolor=color,
            markeredgewidth=1.0,
            linestyle="none",
            zorder=4,
        )
    elif low is None or high is None:
        unknown_ax.plot(
            0.5,
            y,
            marker=marker,
            markersize=size,
            markerfacecolor=color if filled else "white",
            markeredgecolor=color,
            linestyle="none",
            zorder=4,
        )


def save_figure(fig: plt.Figure, out: Path, name: str) -> list[Path]:
    paths = [out / f"{name}.pdf", out / f"{name}.png"]
    fig.savefig(paths[0], bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(paths[1], dpi=180, bbox_inches="tight")
    plt.close(fig)
    return paths


def plot_depth_overview(
    libraries: list[dict[str, str]],
    host_summary: list[dict[str, str]],
    out: Path,
) -> list[Path]:
    by_host = {r["species_current"]: r for r in host_summary}
    hosts = sorted(by_host, key=lambda host: (not host.startswith("Physalia "), host))
    y_by_host = {
        host: i + (1.3 if not host.startswith("Physalia ") else 0)
        for i, host in enumerate(hosts)
    }
    nrows = max(y_by_host.values()) + 1
    fig = plt.figure(figsize=(9.5, 12.5))
    ax = fig.add_axes((0.29, 0.10, 0.32, 0.80))
    unknown_ax = fig.add_axes((0.62, 0.10, 0.065, 0.80), sharey=ax)
    values_ax = fig.add_axes((0.71, 0.10, 0.27, 0.80), sharey=ax)
    fig.text(
        0.025,
        0.97,
        "Collection depth and associate detections across all 205 specimens",
        fontsize=12,
        weight="bold",
    )
    fig.text(
        0.025,
        0.938,
        "Depth marks show individual specimens below the Physalia host-group summaries; fractions show accepted detections.",
        fontsize=9,
    )
    for host in hosts:
        y = y_by_host[host]
        rows = sorted(
            (r for r in libraries if r["species_current"] == host), key=depth_sort
        )
        n = int(by_host[host]["total_libraries"])
        if len(rows) != n:
            raise ValueError(f"Host denominator differs between input tables: {host}")
        if hosts.index(host) % 2 == 0:
            for axis in (ax, unknown_ax, values_ax):
                axis.axhspan(y - 0.47, y + 0.47, color=PALE, zorder=0)
        if host.startswith("Physalia "):
            ax.plot(0, y, marker="^", color=ORANGE, markersize=6, linestyle="none")
        else:
            known = [row for row in rows if row["depth_status"] != "unknown"]
            offsets = np.linspace(-0.31, 0.31, len(known)) if len(known) > 1 else [0]
            for row, offset in zip(known, offsets):
                draw_depth(ax, unknown_ax, row, y + offset, size=3.4)
            missing = len(rows) - len(known)
            if missing:
                unknown_ax.text(
                    0.5, y, str(missing), ha="center", va="center", fontsize=10
                )
        for column, (field, _) in enumerate(DOMAIN_FIELDS):
            present = int(by_host[host][f"{field}_positive_libraries"])
            fraction = present / n
            values_ax.add_patch(
                Rectangle(
                    (column - 0.45, y - 0.42),
                    0.9,
                    0.84,
                    color=plt.cm.Blues(0.05 + 0.7 * fraction),
                    linewidth=0,
                )
            )
            values_ax.text(
                column,
                y,
                f"{present}/{n}",
                ha="center",
                va="center",
                fontsize=9,
                color="white" if fraction >= 0.7 else DARK,
            )
    ax.set_yticks(
        [y_by_host[h] for h in hosts],
        [f"{h}  ({by_host[h]['total_libraries']})" for h in hosts],
        fontsize=9,
        fontstyle="italic",
    )
    ax.set_ylim(nrows - 0.35, -1.6)
    depth_axis(ax, 1600)
    ax.text(
        0, -1.1, "Physalia: all 151 at 0 m", fontsize=10, color=ORANGE, weight="bold"
    )
    ax.text(0, 5.0, "150 assigned surface; 1 recorded", fontsize=9, color=ORANGE)
    ax.text(0, 5.6, "Other siphonophores: 54 specimens", fontsize=9, weight="bold")
    unknown_ax.set_xlim(0, 1)
    unknown_ax.set_xticks([])
    unknown_ax.text(
        0.5, -1.0, "Unknown\ndepth (n)", ha="center", va="center", fontsize=8
    )
    unknown_ax.tick_params(left=False, labelleft=False)
    unknown_ax.spines[:].set_visible(False)
    values_ax.set_xlim(-0.5, 3.5)
    values_ax.set_xticks([])
    values_ax.tick_params(left=False, labelleft=False)
    values_ax.spines[:].set_visible(False)
    for column, (_, label) in enumerate(DOMAIN_FIELDS):
        values_ax.text(column, -1, label, ha="center", fontsize=8, weight="bold")
    fig.text(0.025, 0.078, "Host labels (number of specimens)", fontsize=10)
    fig.legend(
        handles=[
            Line2D(
                [], [], color=BLUE, marker="o", linestyle="none", label="Recorded depth"
            ),
            Line2D(
                [],
                [],
                color=BLUE,
                marker="|",
                linewidth=1.5,
                label="Recorded depth interval",
            ),
            Line2D(
                [],
                [],
                color=ORANGE,
                marker="^",
                linestyle="none",
                label="Physalia surface convention (marks summarize host groups)",
            ),
        ],
        loc="lower center",
        bbox_to_anchor=(0.52, 0.030),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    fig.text(
        0.025,
        0.015,
        "Depth intervals are shown without imputed midpoints. Detection categories use different reference catalogs and are not comparable richness estimates.",
        fontsize=9,
    )
    return save_figure(fig, out, "depth_host_overview")


def plot_dt68(rows: list[dict[str, str]], out: Path) -> list[Path]:
    rows = sorted(rows, key=lambda r: (r["species_current"], depth_sort(r)))
    fig = plt.figure(figsize=(9.5, 8.7))
    ax = fig.add_axes((0.285, 0.16, 0.295, 0.70))
    unknown_ax = fig.add_axes((0.585, 0.16, 0.06, 0.70), sharey=ax)
    context_ax = fig.add_axes((0.655, 0.16, 0.325, 0.70), sharey=ax)
    fig.text(
        0.025,
        0.965,
        "DT-68 detections in Nanomia and Resomia",
        fontsize=15,
        weight="bold",
    )
    fig.text(
        0.025,
        0.925,
        "All 24 specimens, including three with unresolved collection depths",
        fontsize=11,
    )
    last_species = None
    for y, row in enumerate(rows):
        if row["species_current"] != last_species:
            if last_species is not None:
                for axis in (ax, unknown_ax, context_ax):
                    axis.axhline(y - 0.5, color="#AAB7C0", linewidth=0.8)
            last_species = row["species_current"]
        if y % 2 == 0:
            for axis in (ax, unknown_ax, context_ax):
                axis.axhspan(y - 0.48, y + 0.48, color=PALE, zorder=0)
        detected = is_true(row["dt68_detected"])
        grade = row["dt68_grade"]
        color = ORANGE if detected else "#667782"
        marker = "s" if grade == "high_confidence" else "o"
        draw_depth(ax, unknown_ax, row, y, color, detected, marker, size=6)
        context_ax.text(
            0.00,
            y,
            row["locality"] or row["ocean_region"] or "Not recorded",
            va="center",
            fontsize=9,
        )
        context_ax.text(
            0.67, y, row["collection_year"] or "—", ha="center", va="center", fontsize=9
        )
        context_ax.text(
            0.91,
            y,
            row["flowcell_count"] or "?",
            ha="center",
            va="center",
            fontsize=9,
            color=DARK,
        )
    ax.set_yticks(
        range(len(rows)),
        [
            f"{r['species_current'].replace('Nanomia', 'N.').replace('Resomia', 'R.')}  ·  {abbreviated_sample(r)}"
            for r in rows
        ],
        fontsize=9,
    )
    ax.set_ylim(len(rows) - 0.35, -1.6)
    depth_axis(ax, 1000)
    unknown_ax.set_xlim(0, 1)
    unknown_ax.set_xticks([])
    unknown_ax.tick_params(left=False, labelleft=False)
    unknown_ax.spines[:].set_visible(False)
    unknown_ax.text(0.5, -1.0, "Depth\nunknown", ha="center", va="center", fontsize=9)
    context_ax.set_xlim(-0.02, 1)
    context_ax.set_xticks([])
    context_ax.tick_params(left=False, labelleft=False)
    context_ax.spines[:].set_visible(False)
    for x, label, align in (
        (0, "Region / locality", "left"),
        (0.67, "Year", "center"),
        (0.91, "Physical\nflowcells", "center"),
    ):
        context_ax.text(
            x, -1.0, label, ha=align, va="center", fontsize=9, weight="bold"
        )
    fig.legend(
        handles=[
            Line2D(
                [],
                [],
                color=ORANGE,
                marker="s",
                linestyle="none",
                label="DT-68: high-confidence detection",
            ),
            Line2D(
                [],
                [],
                color=ORANGE,
                marker="o",
                linestyle="none",
                label="DT-68: validated detection",
            ),
            Line2D(
                [],
                [],
                markeredgecolor="#667782",
                markerfacecolor="white",
                marker="o",
                linestyle="none",
                label="DT-68 not detected",
            ),
        ],
        loc="lower center",
        bbox_to_anchor=(0.50, 0.075),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    fig.text(
        0.025,
        0.047,
        "Both detection grades satisfy read-mapping criteria; high-confidence detections also satisfy same-library assembly criteria.",
        fontsize=9,
    )
    fig.text(
        0.025,
        0.025,
        "Depth intervals have no imputed midpoint. All libraries are shown; flowcell-based inference requires one physical flowcell.",
        fontsize=9,
    )
    return save_figure(fig, out, "dt68_depth_context")


def plot_tissue_context(
    tissue_rows: list[dict[str, str]],
    incidence_rows: list[dict[str, str]],
    out: Path,
) -> list[Path]:
    """Show sampled tissue denominators and paired regional detection fractions."""
    tissue_rows = sorted(
        tissue_rows,
        key=lambda row: (
            row["host_group"] != "Physalia",
            -int(row["total_libraries"]),
            row["tissue_group"],
        ),
    )
    fig = plt.figure(figsize=(10, 9.8))
    bar_ax = fig.add_axes((0.235, 0.615, 0.325, 0.245))
    values_ax = fig.add_axes((0.63, 0.615, 0.32, 0.245), sharey=bar_ax)
    heat_ax = fig.add_axes((0.235, 0.16, 0.73, 0.30))
    fig.text(
        0.025,
        0.966,
        "Sampled tissues and regional Mycoplasmatales detections",
        fontsize=15,
        weight="bold",
    )
    fig.text(
        0.025,
        0.91,
        "A  Tissue sampling across all 205 specimens",
        fontsize=12,
        weight="bold",
    )
    for i, row in enumerate(tissue_rows):
        count = int(row["total_libraries"])
        color = ORANGE if row["host_group"] == "Physalia" else BLUE
        bar_ax.barh(i, count, color=color, height=0.70)
        bar_ax.text(count + 2, i, str(count), va="center", fontsize=10)
        for j, (field, _) in enumerate(DOMAIN_FIELDS):
            present = int(row[f"{field}_positive_libraries"])
            fraction = present / count
            values_ax.add_patch(
                Rectangle(
                    (j - 0.44, i - 0.42),
                    0.88,
                    0.84,
                    color=plt.cm.Blues(0.05 + 0.7 * fraction),
                    linewidth=0,
                )
            )
            values_ax.text(
                j,
                i,
                f"{present}/{count}",
                ha="center",
                va="center",
                fontsize=10,
                color="white" if fraction >= 0.7 else DARK,
            )
    bar_ax.set_yticks(
        range(len(tissue_rows)),
        [
            f"{'Physalia' if r['host_group'] == 'Physalia' else 'Other hosts'}: {'tissue not recorded' if r['tissue_group'] == 'unknown' else r['tissue_group']}"
            for r in tissue_rows
        ],
        fontsize=10,
    )
    bar_ax.set_ylim(len(tissue_rows) - 0.4, -0.6)
    bar_ax.set_xlim(0, max(int(r["total_libraries"]) for r in tissue_rows) * 1.13)
    bar_ax.set_xlabel("Specimens", fontsize=10)
    bar_ax.tick_params(axis="y", length=0)
    bar_ax.spines["left"].set_visible(False)
    values_ax.set_xlim(-0.5, 3.5)
    values_ax.set_xticks(range(4), [label for _, label in DOMAIN_FIELDS], fontsize=10)
    values_ax.xaxis.tick_top()
    values_ax.tick_params(axis="x", length=0, pad=8)
    values_ax.tick_params(axis="y", left=False, labelleft=False)
    values_ax.spines[:].set_visible(False)
    fig.text(0.63, 0.887, "Specimens with detections / specimens sampled", fontsize=10)
    selected = [
        r
        for r in incidence_rows
        if r["scope"] == "all_libraries" and r["ocean_region"] != "all_regions"
    ]
    regions = sorted({r["ocean_region"] for r in selected})
    lookup = {(r["target"], r["tissue"], r["ocean_region"]): r for r in selected}
    target_rows = [
        (target, tissue)
        for target in ("MAGSP0005", "MAGSP0007", "MAGSP0010")
        for tissue in ("all_tissues", "tentacle")
    ]
    values = np.ma.masked_all((len(target_rows), len(regions)))
    for i, (target, tissue) in enumerate(target_rows):
        for j, region in enumerate(regions):
            row = lookup[(target, tissue, region)]
            denominator = int(row["total_libraries"])
            if denominator:
                values[i, j] = float(row["fraction_positive"])
    cmap = plt.cm.Blues.copy()
    cmap.set_bad("#E0E0E0")
    im = heat_ax.imshow(values, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i, (target, tissue) in enumerate(target_rows):
        for j, region in enumerate(regions):
            row = lookup[(target, tissue, region)]
            denominator = int(row["total_libraries"])
            text = f"{row['positive_libraries']}/{denominator}" if denominator else "—"
            fraction = float(row["fraction_positive"]) if denominator else 0
            heat_ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                fontsize=9,
                color="white" if fraction >= 0.6 else DARK,
            )
    heat_ax.set_yticks(
        range(len(target_rows)),
        [
            f"{target}  ·  {'all tissues' if tissue == 'all_tissues' else 'tentacle'}"
            for target, tissue in target_rows
        ],
        fontsize=10,
    )
    region_labels = {
        "Gulf of California": "Gulf of\nCalifornia",
        "Gulf of Mexico": "Gulf of\nMexico",
        "Central Pacific": "Central\nPacific",
    }
    heat_ax.set_xticks(
        range(len(regions)),
        [region_labels.get(r, r.replace(" ", "\n", 1)) for r in regions],
        fontsize=9,
    )
    heat_ax.tick_params(length=0)
    for y in (1.5, 3.5):
        heat_ax.axhline(y, color="white", linewidth=2)
    heat_ax.set_xticks(np.arange(-0.5, len(regions), 1), minor=True)
    heat_ax.grid(which="minor", color="white", linewidth=0.6)
    heat_ax.tick_params(which="minor", bottom=False)
    fig.text(
        0.025,
        0.54,
        "B  Regional detections among Physalia: all tissues and tentacles only",
        fontsize=12,
        weight="bold",
    )
    fig.text(
        0.235,
        0.49,
        "Fractions give detected specimens / sampled specimens; all 151 Physalia contribute to the all-tissue rows.",
        fontsize=10,
    )
    colorbar_ax = fig.add_axes((0.73, 0.058, 0.23, 0.012))
    fig.colorbar(
        im,
        cax=colorbar_ax,
        orientation="horizontal",
        ticks=[0, 0.5, 1],
        label="Fraction detected",
    )
    fig.text(
        0.025,
        0.07,
        "A dash indicates no sampled specimen in that tissue–region group.",
        fontsize=9,
    )
    fig.text(
        0.025,
        0.043,
        "These are descriptive comparisons; host species and collection conditions also vary among regions.",
        fontsize=9,
    )
    return save_figure(fig, out, "tissue_context")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT.resolve()))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir", type=Path, default=ROOT / "figures/collection_context"
    )
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "figures/collection_context"
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    source_paths = {
        name: args.input_dir / name
        for name in (
            "library_context.tsv",
            "host_detection_summary.tsv",
            "dt68_specimens.tsv",
            "tissue_summary.tsv",
            "tissue_region_incidence.tsv",
        )
    }
    tables = {name: read_tsv(path) for name, path in source_paths.items()}
    if len(tables["library_context.tsv"]) != 205:
        raise ValueError("Review figure scope when the 205-library cohort changes")
    outputs = plot_depth_overview(
        tables["library_context.tsv"],
        tables["host_detection_summary.tsv"],
        args.output_dir,
    )
    outputs += plot_dt68(tables["dt68_specimens.tsv"], args.output_dir)
    outputs += plot_tissue_context(
        tables["tissue_summary.tsv"],
        tables["tissue_region_incidence.tsv"],
        args.output_dir,
    )
    script = Path(__file__).resolve()
    snapshot = args.output_dir / "provenance/plot_collection_context.py"
    snapshot.parent.mkdir(exist_ok=True)
    shutil.copyfile(script, snapshot)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    manifest = {
        "status": "complete",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_revision": revision,
        "working_directory": str(ROOT),
        "source_git_status_porcelain": subprocess.check_output(
            ["git", "status", "--porcelain"], cwd=ROOT, text=True
        ).splitlines(),
        "command": shlex.join([sys.executable, *sys.argv]),
        "source_snapshot": relative(snapshot),
        "script_sha256": sha256(script),
        "inputs_sha256": {
            relative(path): sha256(path) for path in source_paths.values()
        },
        "outputs_sha256": {relative(path): sha256(path) for path in outputs},
        "environment": {
            "python": platform.python_version(),
            "matplotlib": matplotlib.__version__,
            "numpy": np.__version__,
        },
        "scope": "All 205 libraries descriptive; only one-physical-flowcell libraries eligible for flowcell inference.",
    }
    (args.output_dir / "plot_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(f"Created {len(outputs) // 2} PDF/PNG figure pairs in {args.output_dir}")


if __name__ == "__main__":
    main()
