#!/usr/bin/env python3
"""Generate descriptive manuscript values from accepted detections and metadata.

No sequence analyses or statistical tests. See docs/dalila_revision.md.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures/dalila_revision"
INPUTS: dict[str, str] = {}
UNRESOLVED = {
    "Nanomia sp. 1",
    "Apolemia sp. 3 CWD-2005",
    "Cordagalma",
    "Forskalia",
    "Physonectae",
    "Praya",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(relative: str) -> list[dict[str, str]]:
    path = ROOT / relative
    INPUTS[relative] = sha256(path)
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name: str, rows: list[dict]) -> None:
    with (OUT / name).open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metadata = read("figures/cohort/sample_metadata.tsv")
    context = read("figures/collection_context/library_context.tsv")
    read("config/publication_host_names.tsv")
    meta = {r["library_id"]: r for r in metadata}
    if len(meta) != len(metadata) or {r["library_id"] for r in context} != set(meta):
        raise ValueError("Metadata and detection identifiers must match uniquely")
    if len(context) != len(metadata) or any(
        r["species_current"] != meta[r["library_id"]]["species_current"]
        for r in context
    ):
        raise ValueError("Context host labels do not match publication metadata")
    counts = Counter(r["species_current"] for r in metadata)
    # Explicit reviewed unresolved list prevents silently classifying a new label.
    if not UNRESOLVED <= counts.keys():
        raise ValueError("Review named-species scope for changed host labels")
    named = set(counts) - UNRESOLVED
    myco = [r for r in context if int(r["metamycoplasmataceae_detected"])]
    values = {
        "PhysaliaN": sum(
            n for host, n in counts.items() if host.startswith("Physalia ")
        ),
        "NanomiaN": sum(n for host, n in counts.items() if host.startswith("Nanomia ")),
        "NanomiaBijugaN": counts["Nanomia bijuga"],
        "NanomiaBijugaBacterialPositive": sum(
            r["species_current"] == "Nanomia bijuga" and int(r["bacterial_detected"])
            for r in context
        ),
        "OtherN": sum(
            n
            for host, n in counts.items()
            if not host.startswith(("Physalia ", "Nanomia "))
        ),
        "NamedSpecimens": sum(counts[host] for host in named),
        "NamedSpecies": len(named),
        "UnresolvedSpecimens": sum(counts[host] for host in UNRESOLVED),
        "UnresolvedGroups": len(UNRESOLVED),
        "MycoHostSpecies": len({r["species_current"] for r in myco} & named),
        "MycoUnresolvedSpecimens": sum(
            r["species_current"] in UNRESOLVED for r in myco
        ),
    }
    physalia = [r for r in context if r["species_current"].startswith("Physalia ")]
    basins = {
        "Pacific": {
            "Central Pacific",
            "NW Pacific",
            "SW Pacific",
            "SE Pacific",
            "Gulf of California",
        },
        "Indian": {"E Indian", "W Indian"},
        "Atlantic": {"SW Atlantic", "NW Atlantic", "NE Atlantic", "Gulf of Mexico"},
    }
    if {r["ocean_region"] for r in physalia} != set.union(*basins.values()):
        raise ValueError("Review basin mapping for changed regions")
    for basin, regions in basins.items():
        subset = [r for r in physalia if r["ocean_region"] in regions]
        values[basin + "CladeA"] = sum(int(r["clade_a_detected"]) for r in subset)
        values[basin + "Physalia"] = len(subset)
    settings = json.loads((ROOT / "config/reporting_sources.json").read_text())
    INPUTS["config/reporting_sources.json"] = hashlib.sha256(Path("config/reporting_sources.json").read_bytes()).hexdigest()
    samples = read(settings["samples"])
    routes = {r["sample_id"]: r["host_route"] for r in samples}
    mags = [
        r
        for r in read(settings["mag_catalog"])
        if r["catalog_status"] == "retained"
    ]
    assemblies = read(settings["assembly_summary"])
    route_rows = []
    for route in sorted(set(routes.values())):
        retained = [r for r in mags if routes[r["sample_id"]] == route]
        route_rows.append(
            {
                "host_route": route,
                "libraries": sum(value == route for value in routes.values()),
                "assembled_libraries": sum(
                    r["host_route"] == route and r["assembly_status"] == "assembled"
                    for r in assemblies
                ),
                "retained_mags": len(retained),
                "mag_yielding_libraries": len({r["sample_id"] for r in retained}),
            }
        )
    write("assembly_routes.tsv", route_rows)
    values["HostDepletedMAGs"] = sum(routes[r["sample_id"]] != "none" for r in mags)
    values["OtherMAGs"] = sum(routes[r["sample_id"]] == "none" for r in mags)
    (OUT / "values.tex").write_text(
        "% Generated by scripts/summarize_dalila_revision.py; do not edit.\n"
        + "".join(
            "\\newcommand{\\dalila" + key + "}{" + str(value) + "}\n"
            for key, value in values.items()
        )
    )
    write(
        "host_species_summary.tsv",
        [
            {
                "host": host,
                "named_species": host in named,
                "specimens": counts[host],
                "mycoplasmataceae_positive": sum(
                    r["species_current"] == host for r in myco
                ),
            }
            for host in sorted(counts)
        ],
    )
    physalia = [r for r in context if r["species_current"].startswith("Physalia ")]
    write(
        "clade_a_regions.tsv",
        [
            {
                "ocean_region": region,
                "physalia_specimens": sum(
                    r["ocean_region"] == region for r in physalia
                ),
                "clade_a_positive": sum(
                    r["ocean_region"] == region and int(r["clade_a_detected"])
                    for r in physalia
                ),
                "positive_library_ids": ";".join(
                    r["library_id"]
                    for r in physalia
                    if r["ocean_region"] == region and int(r["clade_a_detected"])
                ),
            }
            for region in sorted({r["ocean_region"] for r in physalia})
        ],
    )
    # Source life-stage labels are inconsistent and are retained only in metadata.
    (OUT / "recorded_life_stages.tsv").unlink(missing_ok=True)
    (OUT / "summary.json").write_text(json.dumps(values, indent=2) + "\n")
    outputs = {
        str(p.relative_to(ROOT)): sha256(p)
        for p in OUT.iterdir()
        if p.name != "provenance.json"
    }
    provenance = {
        "status": "complete",
        "source_revision": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_note": "Working-tree reporting changes are identified by script and input checksums.",
        "script_sha256": sha256(Path(__file__)),
        "inputs_sha256": INPUTS,
        "outputs_sha256": outputs,
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(values, indent=2))


if __name__ == "__main__":
    main()
