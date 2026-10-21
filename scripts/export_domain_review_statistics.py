#!/usr/bin/env python3
"""Generate affected manuscript values from the corrected statistical tables."""

import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


def main():
    out = Path("figures/catalog_review")
    inputs = {}

    def read(path):
        path = Path(path)
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    settings = json.loads(Path("config/reporting_sources.json").read_text())
    inputs["config/reporting_sources.json"] = hashlib.sha256(Path("config/reporting_sources.json").read_bytes()).hexdigest()
    root = Path(settings["analysis_root"])
    grades = read(root / "grades/bacterial_grades.tsv")
    real = [r for r in grades if r["role"] != "decoy"]
    present = [r for r in real if r["grade"] in {"validated", "high_confidence"}]
    samples = {r["sample_id"]: r for r in read(settings["samples"])}
    values = {
        "Detections": len(present),
        "AssemblySupported": sum(r["grade"] == "high_confidence" for r in present),
        "TracePairs": sum(r["grade"] == "trace" for r in real),
    }
    for label, depleted in [("Depleted", True), ("Full", False)]:
        rows = [
            r
            for r in present
            if (samples[r["sample_id"]]["host_route"] != "none") == depleted
        ]
        values[label + "Detections"] = len(rows)
        values[label + "Supported"] = sum(r["grade"] == "high_confidence" for r in rows)
    tests = read(root / "grades/contamination_tests.tsv")
    statuses = Counter(r["status"] for r in tests)
    values.update(
        BatchTargets=len(tests),
        BatchTested=statuses["tested"],
        BatchInsufficient=statuses["below_minimum_eligible_presences"],
        BatchConstant=statuses["untestable_constant_statistic"],
        BatchFlagged=sum(r["probable_contaminant"] == "true" for r in tests),
    )
    threshold_rows = []
    for label, threshold, folder, grade in [
        ("Five", 5, "sensitivity/breadth_5pct", "grade_at_5pct"),
        ("Ten", 10, "primary", "grade"),
        ("Twenty", 20, "sensitivity/breadth_20pct", "grade_at_20pct"),
    ]:
        calls = [r for r in real if r[grade] in {"validated", "high_confidence"}]
        region = next(
            r
            for r in read(root / folder / "physalia_permanova.tsv")
            if r["term"] == "ocean_region"
        )
        permutation = [
            r
            for r in read(root / folder / "physalia_region_permutation.tsv")
            if r["status"] == "tested"
        ]
        significant = sum(float(r["q_value"]) < 0.05 for r in permutation)
        n_genomes = len({r["target_id"] for r in calls})
        n_samples = len({r["sample_id"] for r in calls})
        values.update(
            {
                label + "Pairs": len(calls),
                label + "Targets": n_genomes,
                label + "Libraries": n_samples,
                label + "CommunityN": region["libraries"],
                label + "RegionPercent": f"{100 * float(region['R2']):.1f}",
                label + "RegionP": f"{float(region['Pr(>F)']):.4f}",
                label + "TaxaTested": len(permutation),
                label + "TaxaSignificant": significant,
            }
        )
        threshold_rows.append(
            f"{threshold} & {len(calls)} & {n_genomes} & {n_samples} & {region['libraries']} & {float(region['R2']):.3f} & {float(region['Pr(>F)']):.4f} & {significant}/{len(permutation)} \\\\"
        )
    (out / "statistics.tex").write_text(
        "".join(
            f"\\newcommand{{\\review{k}}}{{{v:,}}}\n"
            if isinstance(v, int)
            else f"\\newcommand{{\\review{k}}}{{{v}}}\n"
            for k, v in values.items()
        )
    )
    (out / "threshold_rows.tex").write_text("\n".join(threshold_rows) + "\n")
    header = r"""\begin{tabular}{rrrrrrrr}
\toprule
\shortstack{Breadth\\(\%)} & Detections & Genomes & Specimens &
\shortstack{Community\\$n$} & \shortstack{Region\\$R^2$} &
\shortstack{Region\\$P$} & \shortstack{Taxa with\\$q<0.05$} \\
\midrule
"""
    (out / "threshold_table.tex").write_text(
        header
        + "\n".join(threshold_rows)
        + "\n"
        + r"\bottomrule"
        + "\n"
        + r"\end{tabular}"
        + "\n"
    )
    (out / "statistics.json").write_text(json.dumps(values, indent=2) + "\n")
    inputs[__file__] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (out / "statistics_provenance.json").write_text(
        json.dumps({"inputs_sha256": inputs}, indent=2) + "\n"
    )
    print(json.dumps(values, indent=2))


if __name__ == "__main__":
    main()
