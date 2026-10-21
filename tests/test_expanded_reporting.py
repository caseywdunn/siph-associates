"""Background flags use the locked read threshold and exclude host/human evidence."""

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "expanded_review",
    Path(__file__).resolve().parents[1] / "scripts/compare_expanded_reporting.py",
)
review = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(review)


def test_background_flags_aggregate_bands_without_counting_human_or_host():
    metadata = [
        {"analysis_sample_id": s, "study": "study", "species_current": "host"}
        for s in ["flagged", "clean"]
    ]
    lineages = [
        ("flagged", "Eukaryota;Metazoa;Arthropoda;Insecta", 97, 20),
        ("flagged", "Eukaryota;Metazoa;Arthropoda;Insecta", 99, 17),
        ("clean", "Eukaryota;Metazoa;Arthropoda;Insecta", 97, 36),
        ("clean", "Eukaryota;Metazoa;Arthropoda;Insecta", 95, 999),
        ("clean", "Eukaryota;Metazoa;Chordata;Mammalia;Homo", 99, 999),
        ("clean", "Eukaryota;Metazoa;Cnidaria;Siphonophorae", 99, 999),
    ]
    rows = [
        dict(sample_id=s, lineage=t, min_identity_band=str(i), read_pairs=str(n))
        for s, t, i, n in lineages
    ]
    flags = review.background_flags(metadata, rows, 97, 37)
    assert flags[0]["background_flag"] == "true"
    assert flags[0]["max_control_pairs"] == 37
    assert flags[1]["background_flag"] == "false"
    assert flags[1]["max_control_pairs"] == 36
    assert flags[1]["Mammalia_pairs"] == 0
