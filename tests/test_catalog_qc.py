"""Protect joins and coordinate accounting used by the expanded catalog review."""

import io
import sys
import tarfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "workflow/scripts"))
from screen_nominations import screen_nominations
from summarize_catalog_qc import covered


@pytest.mark.parametrize(
    "intervals, expected",
    [([], 0), ([(0, 10), (5, 15), (7, 8)], 15), ([(20, 25), (0, 10)], 15)],
)
def test_coding_fraction_counts_overlapping_genes_once(intervals, expected):
    assert covered(intervals) == expected


def test_screen_sources_preserve_support_and_library_identity(tmp_path):
    bracken = tmp_path / "bracken.tsv"
    bracken.write_text(
        "name\ttaxonomy_id\tnew_est_reads\tfraction_total_reads\tkraken_assigned_reads\n"
        "Example bacterium\t123\t40\t0.04\t35\n"
    )
    sylph = tmp_path / "sylph.tsv"
    sylph.write_text(
        "Genome_file\tTaxonomic_abundance\tSequence_abundance\tAdjusted_ANI\tTrue_cov\tkmers_reassigned\n"
        "/db/genome.fa\t2\t3\t97\t4\t5\n"
    )
    nested = io.BytesIO()
    with tarfile.open(fileobj=nested, mode="w:gz") as archive:
        # A comma in the lineage must not be confused with the final count.
        content = b"Bacteria;Example, lineage,7\n"
        member = tarfile.TarInfo("sample.phyloFlash.NTUfull_abundance.csv")
        member.size = len(content)
        archive.addfile(member, io.BytesIO(content))
    phylo = tmp_path / "phylo.tar.gz"
    with tarfile.open(phylo, "w:gz") as archive:
        content = nested.getvalue()
        member = tarfile.TarInfo("sample.phyloFlash.tar.gz")
        member.size = len(content)
        archive.addfile(member, io.BytesIO(content))
    fields, rows = screen_nominations(
        ["library"], {"library": "study"}, [bracken], [bracken], [sylph], [phylo]
    )
    assert len(rows) == 4
    assert all(r["sample_id"] == "library" and r["study"] == "study" for r in rows)
    assert all(set(r) <= set(fields) for r in rows)
    assert [(r["source"], r["rank"], r["raw_support"]) for r in rows] == [
        ("Bracken", "genus", "35"),
        ("Bracken", "species", "35"),
        ("sylph", "genome", "5"),
        ("phyloFlash", "SSU_lineage", "7"),
    ]
    assert rows[0]["estimated_reads"] == "40"
    assert rows[2]["adjusted_ani"] == "97"
    assert rows[3]["candidate_name"] == "Example, lineage"


def test_missing_library_input_is_not_silently_omitted():
    with pytest.raises(ValueError, match="zip"):
        screen_nominations(["library"], {"library": "study"}, [], [], [], [])
