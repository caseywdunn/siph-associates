"""Regression checks for physical-flowcell inference and threshold selection."""

import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = ROOT / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
from phase6_metadata import normalize_manifest, physical_flowcells

spec = importlib.util.spec_from_file_location("contamination_script", SCRIPTS / "test_contamination_flowcell.py")
contamination_script = importlib.util.module_from_spec(spec)
spec.loader.exec_module(contamination_script)


def write_table(path, records, delimiter="\t"):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(records)


def read_table(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_lanes_collapse_but_distinct_physical_flowcells_do_not():
    assert physical_flowcells("I:1:F:L001;I:1:F:L002") == ("I:1:F",)
    assert physical_flowcells("I:1:F:L002;I:2:G:L001") == ("I:1:F", "I:2:G")
    with pytest.raises(ValueError):
        physical_flowcells("I/1/F/L001")


def test_normalization_retains_multi_flowcell_libraries(tmp_path):
    path = tmp_path / "manifest.csv"
    common = {"include_primary": "true", "study": "Study", "species_current": "Physalia physalis",
              "ocean_region": "Atlantic"}
    write_table(path, [
        {**common, "library_id": "Study:A", "sequencing_batches": "I:1:F:L001;I:1:F:L002"},
        {**common, "library_id": "Study:B", "sequencing_batches": "I:1:F:L001;I:2:G:L001"},
    ], ",")
    rows = normalize_manifest(path)
    assert len(rows) == 2
    assert rows[0]["inference_eligible"] == "true"
    assert rows[0]["flowcell"] == "I:1:F"
    assert rows[1]["inference_eligible"] == "false"
    assert rows[1]["inference_exclusion"] == "multiple_flowcells"
    assert rows[1]["flowcell_count"] == 2


def test_contamination_uses_eligible_libraries_and_selected_threshold():
    metadata = [{"sample_id": str(i), "host_species": "Nanomia septata", "ocean_region": "Pacific",
                 "inference_eligible": "true" if i < 8 else "false", "flowcell": f"I:1:{i // 4}"}
                for i in range(10)]
    grades = [{"sample_id": str(i), "target_id": "BAC001", "role": "reference", "taxonomy": "g__A",
               "grade": "validated" if i in (0, 1, 2, 8, 9) else "none",
               "grade_at_20pct": "validated" if i in (0, 8, 9) else "none"} for i in range(10)]
    rule = {"min_flowcells_in_stratum": 2, "min_validated_presences": 3, "permutations": 99,
            "seed": 123, "fdr": 0.05, "identity_annotation_genera": []}
    result = contamination_script.contamination_results(rule, grades, metadata)[0]
    assert result["validated_presences"] == 5
    assert result["eligible_presences"] == result["permutable_presences"] == 3
    assert result["excluded_presences"] == 2
    assert result["status"] == "tested"
    # Neither input record order nor hash-dependent set ordering changes the random stream.
    assert contamination_script.contamination_results(rule, grades[::-1], metadata[::-1])[0] == result
    strict = contamination_script.contamination_results(rule, grades, metadata, "grade_at_20pct")[0]
    assert strict["validated_presences"] == 3
    assert strict["eligible_presences"] == 1
    assert strict["status"] == "below_minimum_eligible_presences"
    assert strict["p_value"] == ""


def test_r_analysis_excludes_multiflowcell_and_uses_threshold(tmp_path):
    config = json.loads((ROOT / "config" / "phase6_analysis.json").read_text())
    rscript = Path(config["software"]["stats_prefix"]) / "bin" / "Rscript"
    if not rscript.exists():
        pytest.skip("configured R environment unavailable")
    config["physalia_models"]["min_presences"] = 3
    config["physalia_models"]["permutations"] = 19
    config_path = tmp_path / "analysis.json"
    config_path.write_text(json.dumps(config))
    metadata, libraries, grades = [], [], []
    for i in range(18):
        sample = f"sample{i:02d}"
        metadata.append({"sample_id": sample, "host_species": "Physalia physalis",
                         "ocean_region": "Atlantic" if i % 2 else "Pacific",
                         "flowcell": f"I:1:F{i // 4}", "inference_eligible": "true" if i < 16 else "false"})
        libraries.append({"sample_id": sample, "input_pairs": 1000 + i * 73})
        for target, hits in (("BAC001", {0, 2, 3, 6, 9, 11, 13, 16, 17}),
                             ("BAC002", {1, 4, 5, 7, 8, 10, 12, 14, 15})):
            grades.append({"sample_id": sample, "target_id": target, "role": "reference",
                           "grade": "validated" if i in hits else "none",
                           "grade_at_5pct": "validated" if i in hits or i == 1 else "none",
                           "grade_at_20pct": "validated" if i in hits and i != 2 else "none"})
    for name, records in (("metadata", metadata), ("libraries", libraries), ("grades", grades)):
        write_table(tmp_path / f"{name}.tsv", records)
    observed = {}
    for column in ("grade", "grade_at_5pct", "grade_at_20pct"):
        outputs = [tmp_path / f"{column}_{name}.tsv" for name in ("design", "models", "permanova", "permutation")]
        subprocess.run([str(rscript), str(SCRIPTS / "fit_physalia_models.R"), str(config_path),
                        str(tmp_path / "grades.tsv"), str(tmp_path / "libraries.tsv"),
                        str(tmp_path / "metadata.tsv"), *map(str, outputs), column, "1"],
                       check=True, capture_output=True, text=True)
        assert sum(int(row["N"]) for row in read_table(outputs[0])) == 16
        permutation = {row["target_id"]: row for row in read_table(outputs[3])}
        assert all(int(row["libraries"]) == 16 for row in permutation.values())
        assert all(row["grade_column"] == column for row in permutation.values())
        observed[column] = int(permutation["BAC001"]["presences"])
    assert observed == {"grade": 7, "grade_at_5pct": 8, "grade_at_20pct": 6}
