"""Reporting annotations for eukaryote evidence without changing classifications."""

from collections import defaultdict


def annotate_detections(records: list[dict]) -> list[dict]:
    """Retain ambiguous ancestors without counting them as extra detections.

    A validated ancestor and validated descendant in one library can represent
    the same organism. Their read assignments remain separate: ambiguous reads
    are never reassigned to a descendant. Trace descendants cannot establish
    presence and therefore cannot suppress a validated ancestor.
    """
    present = defaultdict(set)
    for record in records:
        if record["grade"] in ("validated", "high_confidence"):
            present[record["sample_id"]].add(record["reporting_unit"])

    annotated = []
    for record in records:
        descendants = sorted(
            unit
            for unit in present[record["sample_id"]]
            if unit.startswith(record["reporting_unit"] + ";")
        )
        validated = record["grade"] in ("validated", "high_confidence")
        status = (
            "unresolved_ancestor"
            if descendants
            else "detection"
            if validated
            else "trace_evidence"
        )
        annotated.append(
            {
                **record,
                "reporting_status": status,
                "count_as_detection": str(validated and not descendants).lower(),
                "overlapping_descendant_units": "|".join(descendants),
            }
        )
    return annotated
