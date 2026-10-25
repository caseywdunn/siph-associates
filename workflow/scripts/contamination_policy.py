"""Apply the author-defined two-tier exclusions without removing evidence rows."""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
POLICY_PATHS = ["config/excluded_libraries.tsv", "config/contaminant_taxa.tsv"]

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

class ContaminationPolicy:
    """Match explicit library IDs, normalized GTDB genera, and reviewed sequences."""
    def __init__(self, root: Path = ROOT, catalog: list[dict[str, str]] | None = None):
        self.excluded = {r["sample_id"]: r for r in read_tsv(root / POLICY_PATHS[0])}
        self.taxa = read_tsv(root / POLICY_PATHS[1])
        self.catalog = {r["catalog_id"]: set(filter(None, [r.get("genome_id", ""), *r.get("cluster_members", "").split(",")])) for r in (catalog or []) if r["catalog_id"]}

    def annotate(self, row: dict[str, str]) -> dict[str, str]:
        sample = row.get("analysis_sample_id") or row.get("sample_id", "")
        if "__" not in sample and row.get("library_id"):
            sample = row["library_id"].replace(":", "__", 1)
        exclusion = self.excluded.get(sample)
        genus = next((x[3:].split("_")[0] for x in row.get("taxonomy", "").split(";") if x.startswith("g__")), "")
        clusters = self.catalog.get(row.get("target_id", ""), {row.get("species_cluster", "")})
        match = None
        for taxon in self.taxa:
            kind = taxon["match_type"]
            if kind == "genus" and genus == taxon["name"]:
                match = taxon
            elif kind == "genome_cluster" and taxon["name"] in clusters:
                match = taxon
            elif kind == "eukaryote_ssu" and sample in taxon["sample_ids"].split(";") and (
                taxon["name"] in row.get("ncbi_top_hit", "") or row.get("reporting_unit") == "Nematozoa"
            ):
                match = taxon
        return {**row,
            "excluded_library": str(exclusion is not None).lower(),
            "sample_exclusion_reason": exclusion["reason"] if exclusion else "",
            "contaminant_taxon": str(match is not None).lower(),
            "contaminant_group": match["category"] if match else "",
            "contaminant_basis": match["reason"] if match else "",
            "biological_include": str(exclusion is None and match is None).lower(),
        }

def retained(row: dict[str, str]) -> bool:
    """Unannotated historical inputs retain their original interpretation."""
    return row.get("biological_include", "true") == "true"
