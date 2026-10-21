"""Current processing routes must not inherit obsolete reference-free labels."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from reporting_manifest import build_reporting_manifest  # noqa: E402


def test_current_route_overrides_baseline_and_preserves_read_identity():
    original = [
        dict(
            library_id="study:CWD12",
            specimen_id="voucher",
            host_reference="none",
            r1_paths="original_R1.fastq.gz",
            species_current="old name",
        )
    ]
    reviewed = [
        dict(
            library_id="study:CWD12",
            analysis_sample_id="study__CWD12",
            specimen_id="edited voucher",
            host_reference="incorrect metadata route",
            r1_paths="incorrect_reads.fastq.gz",
            species_current="Halistemma rubrum",
        )
    ]
    cohort = [dict(sample_id="study__CWD12", host_route="H_rubrum")]
    result = build_reporting_manifest(original, reviewed, cohort)
    assert result == [
        dict(
            original[0], host_reference="H_rubrum", species_current="Halistemma rubrum"
        )
    ]
    assert original[0]["host_reference"] == "none"


@pytest.mark.parametrize("cohort_ids", [["other"], ["sample", "sample"]])
def test_route_join_rejects_mismatched_or_duplicate_cohort(cohort_ids):
    with pytest.raises(ValueError):
        build_reporting_manifest(
            [dict(library_id="library", host_reference="none")],
            [dict(library_id="library", analysis_sample_id="sample")],
            [dict(sample_id=s, host_route="H_rubrum") for s in cohort_ids],
        )
