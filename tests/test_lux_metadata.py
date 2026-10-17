"""Museum comparisons preserve distinctions between measurements and context."""

import sys
from decimal import Decimal
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from review_lux_metadata import coordinate_status, depth_interval, museum_date


def test_depth_requires_explicit_units_and_preserves_zero():
    assert depth_interval("0-50 meters") == (Decimal(0), Decimal(50))
    assert depth_interval("387-387 meters") == (Decimal(387), Decimal(387))
    assert depth_interval("1289.68 meters") == (Decimal("1289.68"),) * 2
    assert depth_interval("2854.16") is None
    assert depth_interval("100 feet") is None
    with pytest.raises(ValueError):
        depth_interval("50-0 meters")


def test_coordinate_comparison_detects_sign_errors_without_imputing_points():
    row = {"longitude": "-71.67", "latitude": "41.34"}
    west = "POLYGON((-71.7 41.3, -71.6 41.3, -71.6 41.4, -71.7 41.3))"
    east = "POLYGON((71.6 41.3, 71.7 41.3, 71.7 41.4, 71.6 41.3))"
    assert coordinate_status(west, row) == "within_locality_envelope"
    assert coordinate_status(east, row) == "outside_locality_envelope"
    assert coordinate_status("", row) == "not_comparable"
    assert coordinate_status(west, {**row, "longitude": ""}) == "not_comparable"
    crossing = "POLYGON((-179 41.3, 179 41.3, 179 41.4, -179 41.3))"
    assert coordinate_status(crossing, row) == "not_comparable"
    assert row == {"longitude": "-71.67", "latitude": "41.34"}


def test_date_extraction_requires_one_unambiguous_explicit_date():
    assert museum_date("Collected 03 Jul 2002 LatLon 40 517 -68 217") == "2002-07-03"
    assert museum_date("ROV Doc Ricketts D1398-SS8") == ""
    assert museum_date("Collected 03 Jul 2002 accessioned 04 Aug 2003") == ""
