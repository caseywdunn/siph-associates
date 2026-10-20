#!/usr/bin/env python3
"""Apply Church's dated responses to the publication overlay and review worksheet."""

import hashlib
import json
from collections import Counter, defaultdict
from decimal import Decimal

from export_sample_metadata import (
    INPUTS,
    PUBLIC_FIELDS,
    ROOT,
    curate_metadata,
    load_metadata,
    read_table,
    write_table,
)

REVIEW = ROOT / "data/metadata/manual-review/church2025-2026-10-19"
SOURCE = REVIEW / "2026-10-19_church2025_review_responses.tsv"
FIELDS = {
    "museum_voucher": "voucher_id",
    "biosample_accession": "biosample",
    "field_event_id": "collection_id",
    "tissue_sample_label": "sample_id",
    "tissue": "tissue",
    "collection_date": "collection_date",
    "coordinates": "coordinates",
}


def main() -> None:
    responses = read_table(SOURCE)
    worksheet_path = (
        ROOT / "data/metadata/manual-review/2026-10-17_church2025_review.tsv"
    )
    worksheet = read_table(worksheet_path)
    overlay_path = ROOT / INPUTS[4]
    overlay = read_table(overlay_path)
    before = {r["library_id"]: r for r in load_metadata()}
    baseline = {
        r["library_id"]: r
        for r in curate_metadata(
            read_table(ROOT / INPUTS[0]), read_table(ROOT / INPUTS[1]), [], []
        )
    }
    answers = {}
    for r in responses:
        key = ("Church2025:" + r["specimen_id"], FIELDS[r["field"]])
        if key in answers or key[0] not in before:
            raise ValueError(f"Duplicate or unmatched author response: {key}")
        answers[key] = r
    requested = {(r["library_id"], r["field"]) for r in worksheet}
    returned = {k for k, v in answers.items() if v["question"] != "extra"}
    if len(responses) != 92 or len(requested) != 70 or requested != returned:
        raise ValueError(
            "Responses must cover the exact 70 requested entries plus 22 corrections"
        )
    patches, note_parts, locators = {}, defaultdict(list), defaultdict(list)
    for (library, field), r in answers.items():
        decision = r["decision"]
        if decision not in {"set_value", "keep_current", "not_assigned", "unresolved"}:
            raise ValueError(f"Unrecognized review decision: {decision}")
        coordinates = field == "coordinates"
        fields = ["latitude", "longitude"] if coordinates else [field]
        current = (
            [v.strip() for v in r["current_value"].split(",")]
            if coordinates
            else [r["current_value"]]
        )
        values = (
            [v.strip() for v in r["reviewed_value"].split(",")]
            if coordinates
            else [r["reviewed_value"]]
        )
        if len(fields) != len(values) or len(fields) != len(current):
            raise ValueError(f"Malformed coordinate pair: {library}")
        locator = f"specimen_id={r['specimen_id']};field={r['field']};question={r['question']}"
        locators[library].append(locator)
        for name, old, value in zip(fields, current, values):
            observed = before[library].get(name, "")
            agrees = (
                Decimal(observed) == Decimal(old) if coordinates else observed == old
            )
            if not agrees:
                raise ValueError(f"Response current value differs: {library} {name}")
            if decision == "keep_current":
                if coordinates and Decimal(value) != Decimal(observed):
                    raise ValueError(
                        f"Keep-current response changes coordinates: {library}"
                    )
                if not coordinates and value != observed:
                    raise ValueError(
                        f"Keep-current response changes field: {library} {name}"
                    )
                value = observed  # Preserve precision/spelling for unchanged values.
            if decision in {"not_assigned", "unresolved"} and (value or observed):
                raise ValueError(
                    f"An unresolved response would erase evidence: {library} {name}"
                )
            patches[(library, name)] = {
                "library_id": library,
                "field": name,
                "original_value": ""
                if name == "notes"
                else baseline[library].get(name, ""),
                "value": value,
                "source": str(SOURCE.relative_to(ROOT)),
                "source_locator": locator,
                "notes": f"Reviewer decision: {decision}. {r['evidence']}. {r['review_notes']}",
            }
        if r["question"] == "extra":
            note = "Submitted sample label uses the institution-prefixed museum identifier; it is not a separate tube number."
        elif decision == "not_assigned":
            note = (
                "No separate tissue/sample label was assigned. "
                + r["review_notes"].split(" No separate tube/sample label")[0]
            )
        elif r["specimen_id"].startswith("OA-"):
            note = "Tissue, museum-voucher deposition and BioSample registration remain unknown; material was provided by Ozren and Alex."
            if r["specimen_id"] == "OA-Bluebottle":
                note += " GON10891 is retained as the collecting-event label, but its use as an event or tissue identifier remains unconfirmed; the sample label is unknown."
        else:
            note = (
                r["review_notes"]
                .replace("Current value", "Previous metadata value")
                .replace("current latitude", "previous metadata latitude")
                .replace("current value.", "previous metadata value.")
            )
        if note not in note_parts[library]:
            note_parts[library].append(note)
    for library, parts in note_parts.items():
        # Existing notes on these records were the unresolved comparisons now reviewed.
        old = {
            r["value"]
            for r in overlay
            if r["library_id"] == library and r["field"] == "notes"
        }
        if any("pending review" not in text for text in old):
            raise ValueError(
                f"Unexpected existing notes require reconciliation: {library}"
            )
        patches[(library, "notes")] = {
            "library_id": library,
            "field": "notes",
            "original_value": "",
            "value": "Author review (Samuel Church, 2026-10-19): " + " ".join(parts),
            "source": str(SOURCE.relative_to(ROOT)),
            "source_locator": " | ".join(locators[library]),
            "notes": "Replaces pending-review wording with the returned decisions; full response and superseded overlay evidence retained.",
        }
    superseded = [r for r in overlay if (r["library_id"], r["field"]) in patches]
    retained = [r for r in overlay if (r["library_id"], r["field"]) not in patches]
    updated = retained + [patches[k] for k in sorted(patches)]
    # Validate the entire curated result before changing the overlay.
    after = curate_metadata(
        read_table(ROOT / INPUTS[0]),
        read_table(ROOT / INPUTS[1]),
        read_table(ROOT / INPUTS[2]),
        updated,
    )
    for row in after:
        for key in (
            "library_id",
            "analysis_sample_id",
            "specimen_id",
            "sra_run",
            "read_pairs",
        ):
            if row[key] != before[row["library_id"]][key]:
                raise ValueError(
                    f"Protected analysis identity changed: {row['library_id']} {key}"
                )
    for row in worksheet:
        response = answers[(row["library_id"], row["field"])]
        for name in ("decision", "reviewed_value", "evidence", "review_notes"):
            row[name] = response[name]
    write_table(
        REVIEW / "metadata_before.tsv",
        list(before.values()),
        PUBLIC_FIELDS,
    )
    write_table(
        REVIEW / "superseded_publication_updates.tsv", superseded, list(overlay[0])
    )
    write_table(overlay_path, updated, list(overlay[0]))
    write_table(worksheet_path, worksheet, list(worksheet[0]))
    summary = {
        "response_entries": len(responses),
        "original_review_entries": len(worksheet),
        "decisions": dict(Counter(r["decision"] for r in responses)),
        "affected_libraries": len(note_parts),
        "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "analysis_inputs_changed": False,
    }
    (REVIEW / "application_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
