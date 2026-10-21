#!/usr/bin/env python3
"""Validate the mapping hold and freeze catalog identities and review accounting."""

import json
import re
import shutil
from collections import defaultdict
from pathlib import Path

from assembly_completion import sha256, table, write_table


def groups(rows: list[dict], kind: str) -> dict[str, set[str]]:
    result = defaultdict(set)
    for row in rows:
        if kind == "mag" and row.get("catalog_status") == "retained":
            result[row["species_cluster"]].add(row["mag_id"])
        elif kind == "bacterial" and row.get("catalog_id"):
            result[row["catalog_id"]].update(row["cluster_members"].split(","))
        elif kind == "viral":
            result[row["votu_id"]].add(re.sub(r"\|qc_\d+_\d+$", "", row["contig_id"]))
    return result


def main() -> None:
    c = json.loads(Path("config/catalog_freeze.json").read_text())
    root = Path(c["work_root"])
    p4 = root / "phase4_catalog/v3"
    hold = json.loads(Path("config/post_expansion_remapping.json").read_text())
    if hold["status"] != "hold" or hold["automatic_remapping"]:
        raise ValueError("Catalog preparation must preserve the explicit mapping hold")
    summary = dict(
        line.split("\t", 1)
        for line in (root / "stages/phase4_catalog_v3.done").read_text().splitlines()
    )
    if summary["status"] != "PASS":
        raise ValueError("Catalog validation has not passed")
    mags = table(root / "phase3_catalog/mags/mag_catalog.tsv")
    decisions = {r["mag_id"]: r for r in table("figures/catalog_qc/bin_decisions.tsv")}
    if any(
        decisions[r["mag_id"]]["decision"] == "exclude"
        for r in mags
        if r["catalog_status"] == "retained"
    ):
        raise ValueError("Excluded domain candidate entered the catalog")
    if any(r["species_cluster"] == "MAGSP0026" for r in mags):
        raise ValueError("Retired false-archaeal identifier was reused")
    out = Path("figures/catalog_qc/catalog_v3")
    out.mkdir(parents=True, exist_ok=True)
    artifacts = [
        p4 / f"{name}{suffix}"
        for name in ("bacterial_catalog.fna", "viral_catalog.fna")
        for suffix in ("", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
    ]
    artifacts += [
        p4 / "bacterial_catalog.manifest.tsv",
        p4 / "viral_associate.manifest.tsv",
        p4 / "viral_endogenous_candidate.manifest.tsv",
        p4 / "crispr/host_links.tsv",
    ]
    artifacts += [
        root / "phase3_catalog/mags/mag_catalog.tsv",
        Path("figures/catalog_qc/bin_decisions.tsv"),
        Path("figures/catalog_qc/contig_decisions.tsv"),
        Path("figures/catalog_qc/curated_mag_quality.tsv"),
        Path("figures/catalog_qc/viral_region_decisions.tsv"),
    ]
    inventory = [
        {"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)}
        for p in artifacts
    ]
    write_table(out / "inventory.tsv", inventory, ["path", "bytes", "sha256"])
    for p in artifacts:
        if p.suffix == ".tsv" and (
            p.parent == p4 or p.name in {"mag_catalog.tsv", "host_links.tsv"}
        ):
            shutil.copyfile(p, out / p.name)
    pairs = [
        (
            "mag",
            c["previous_mag_catalog"],
            root / "phase3_catalog/mags/mag_catalog.tsv",
        ),
        (
            "bacterial",
            c["previous_bacterial_manifest"],
            p4 / "bacterial_catalog.manifest.tsv",
        ),
    ]
    pairs += [
        ("viral", old, root / f"phase3_catalog/votus/{kind}.votus.tsv")
        for kind, old in c["previous_votu_catalogs"].items()
    ]
    crosswalk = []
    for kind, old, new in pairs:
        previous, current = groups(table(old), kind), groups(table(new), kind)
        for old_id, members in previous.items():
            matches = [
                (key, len(members & value))
                for key, value in current.items()
                if members & value
            ]
            for new_id, shared in matches or [("", 0)]:
                crosswalk.append(
                    {
                        "kind": kind,
                        "previous_id": old_id,
                        "current_id": new_id,
                        "shared_members": shared,
                        "status": "retained_id"
                        if new_id == old_id
                        else "changed_cluster"
                        if new_id
                        else "not_retained",
                    }
                )
        old_members = set().union(*previous.values())
        for key, value in current.items():
            if not value & old_members:
                crosswalk.append(
                    {
                        "kind": kind,
                        "previous_id": "",
                        "current_id": key,
                        "shared_members": 0,
                        "status": "new",
                    }
                )
    write_table(out / "catalog_id_crosswalk.tsv", crosswalk, list(crosswalk[0]))
    summary.update(
        mapping_status="held_not_started",
        catalog_root=str(p4),
        curated_mags=len([r for r in mags if r["catalog_status"] == "retained"]),
        mag_species=len(
            {r["species_cluster"] for r in mags if r["catalog_status"] == "retained"}
        ),
        inventory_sha256=sha256(out / "inventory.tsv"),
    )
    (out / "catalog_ready.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
