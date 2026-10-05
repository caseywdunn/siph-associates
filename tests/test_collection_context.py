"""Scientific join, scope, interval and counting checks for collection context."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "collection_context",
    Path(__file__).resolve().parents[1] / "scripts/analyze_collection_context.py",
)
context = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(context)


def specimen(library="Study:A", species="Nanomia septata", **updates):
    row = {
        "library_id": library,
        "species_current": species,
        "depth_m": "10",
        "depth_min_m": "",
        "depth_max_m": "",
        "collection_depth_basis": "source_record",
        "tissue": "siphosome",
        "ocean_region": "Pacific",
        "study": "Study",
    }
    row.update(updates)
    return row


def inputs(rows):
    return {
        "manifest": rows,
        "flowcells": [
            {
                "sample_id": row["library_id"].replace(":", "__"),
                "flowcells": "x;y",
                "flowcell_count": "2",
                "inference_eligible": "false",
            }
            for row in rows
        ],
        "bacteria": [],
        "viruses": [],
        "eukaryotes": [],
        "lineages": [
            {"catalog_id": catalog, "species_cluster": genome, "lineage": lineage}
            for catalog, genome, lineage in [
                ("BAC1", "MAGSP0005", "Clade A"),
                ("BAC2", "MAGSP0007", "Clade B"),
                ("BAC3", "MAGSP0010", "Clade B"),
                ("BAC4", "MAGSP0011", "DT-68"),
            ]
        ],
    }


def test_duplicate_or_ambiguous_canonical_identifiers_fail():
    with pytest.raises(ValueError, match="duplicate"):
        context.index_manifest([specimen(), specimen()])
    with pytest.raises(ValueError, match="Ambiguous"):
        context.index_manifest([specimen("Study:A__B")])


def test_unmatched_evidence_does_not_silently_disappear():
    data = inputs([specimen()])
    data["viruses"] = [{"sample_id": "Study__Missing", "target_id": "Virus1"}]
    with pytest.raises(ValueError, match="Unmatched"):
        context.derive_context(data)


def test_exact_duplicates_count_once_and_disagreement_fails():
    meta = context.index_manifest([specimen()])
    row = {"sample_id": "Study__A", "target_id": "x", "grade": "validated"}
    assert len(context.unique_evidence([row, dict(row)], "target_id", meta)) == 1
    with pytest.raises(ValueError, match="Conflicting"):
        context.unique_evidence(
            [row, {**row, "grade": "high_confidence"}], "target_id", meta
        )


def test_nondetections_and_multiflowcell_specimens_remain_in_denominator():
    data = inputs([specimen(), specimen("Study:B")])
    result, evidence = context.derive_context(data)
    assert len(result) == 2
    assert evidence == []
    assert all(
        row["bacterial_targets"] == row["eukaryote_groups"] == 0 for row in result
    )
    assert context.summarize(result, ())[0]["total_libraries"] == 2
    assert context.summarize(result, ())[0]["single_flowcell_libraries"] == 0


def test_depth_range_stays_interval_and_boundary_crossing_is_separate():
    row = specimen(depth_m="", depth_min_m="0", depth_max_m="20")
    derived = context.depth_context(row)
    assert derived == {
        "depth_status": "recorded_range",
        "depth_group": "non_Physalia_0_20m",
    }
    assert row["depth_m"] == ""
    assert (
        context.depth_context({**row, "depth_max_m": "30"})["depth_group"]
        == "non_Physalia_crosses20m"
    )
    assert (
        context.depth_context({**row, "depth_m": "10"})["depth_status"]
        == "recorded_point_and_range"
    )
    with pytest.raises(ValueError, match="Inconsistent"):
        context.depth_context({**row, "depth_m": "30"})


def test_surface_assignment_is_explicit_and_not_a_nonphysalia_rule():
    row = specimen(
        species="Physalia utriculus",
        depth_m="0",
        collection_depth_basis="curator_assigned_surface",
    )
    assert context.depth_context(row)["depth_status"] == "assigned_surface"
    with pytest.raises(ValueError, match="Surface assignment"):
        context.depth_context({**row, "species_current": "Nanomia septata"})


def test_vibrio_genus_is_not_all_vibrionaceae_and_archaea_stay_separate():
    data = inputs([specimen()])
    data["bacteria"] = [
        {
            "sample_id": "Study__A",
            "target_id": target,
            "grade": "validated",
            "role": "reference",
            "taxonomy": taxonomy,
        }
        for target, taxonomy in [
            ("a", "d__Archaea;p__Thermoproteota"),
            ("b", "d__Bacteria;f__Vibrionaceae;g__Photobacterium"),
            ("c", "d__Bacteria;f__Vibrionaceae;g__Vibrio"),
        ]
    ]
    row = context.derive_context(data)[0][0]
    assert (
        row["bacterial_targets"],
        row["archaeal_targets"],
        row["vibrio_targets"],
    ) == (2, 1, 1)


def test_eukaryote_reporting_flags_and_nonexclusive_roles():
    data = inputs([specimen()])
    data["eukaryotes"] = [
        {
            "sample_id": "Study__A",
            "reporting_unit": unit,
            "grade": grade,
            "role": role,
            "count_as_detection": flag,
        }
        for unit, grade, role, flag in [
            ("Chordata", "validated", "prey", "true"),
            ("Trematoda", "high_confidence", "parasite", "true"),
            ("Trematoda ancestor", "validated", "parasite", "false"),
            ("Trace", "trace", "prey", "false"),
        ]
    ]
    row = context.derive_context(data)[0][0]
    assert row["eukaryote_groups"] == 2
    assert row["prey_detected"] == row["parasite_detected"] == 1
    data["eukaryotes"][-1]["count_as_detection"] = "true"
    with pytest.raises(ValueError, match="conflicts with grade"):
        context.derive_context(data)


def test_tentacle_regional_denominators_include_negative_libraries():
    data = inputs(
        [
            specimen("Study:A", "Physalia utriculus", depth_m="0", tissue="tentacle"),
            specimen("Study:B", "Physalia utriculus", depth_m="0", tissue="tentacle"),
            specimen("Study:C", "Physalia utriculus", depth_m="0", tissue=""),
        ]
    )
    data["bacteria"] = [
        {
            "sample_id": "Study__A",
            "target_id": "BAC1",
            "grade": "validated",
            "role": "mag",
            "taxonomy": "d__Bacteria;f__Metamycoplasmataceae",
        }
    ]
    rows, _ = context.derive_context(data)
    incidence = context.regional_incidence(rows)
    counts = {
        row["tissue"]: (row["positive_libraries"], row["total_libraries"])
        for row in incidence
        if row["scope"] == "all_libraries"
        and row["ocean_region"] == "all_regions"
        and row["target"] == "MAGSP0005"
    }
    assert counts == {"all_tissues": (1, 3), "tentacle": (1, 2)}
