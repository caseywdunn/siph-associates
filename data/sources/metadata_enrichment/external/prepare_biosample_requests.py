"""Join public run accessions to BioSamples and prepare explicit fetch URLs.

Run from the analysis repository root. Headerless RunInfo files are read as
CSV records; only exact run and BioSample accession cells supply the join.
This script does not perform network requests or modify the manifest.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


BASE = Path(__file__).resolve().parent


def main() -> None:
    """Prepare a validated, inspectable accession mapping and batch requests."""
    source_files = sorted(Path("data/sources").glob("*runinfo*"))
    by_run: dict[str, str] = {}
    sources: dict[str, set[str]] = {}
    for path in source_files:
        for row in csv.reader(path.open()):
            if not row or not re.fullmatch(r"[SED]RR\d+", row[0]):
                continue
            samples = {cell for cell in row if re.fullmatch(r"SAM[NED][A-Z]?\d+", cell)}
            if len(samples) != 1:
                raise ValueError(f"Ambiguous BioSample cells: {path}, {row[0]}")
            sample = next(iter(samples))
            if row[0] in by_run and by_run[row[0]] != sample:
                raise ValueError(f"Inconsistent BioSample for run {row[0]}")
            by_run[row[0]] = sample
            sources.setdefault(row[0], set()).add(str(path))

    rows = []
    for library in csv.DictReader(Path("manifest.csv").open()):
        runs = re.findall(r"[SED]RR\d+", library["sra_run"])
        samples = {by_run[run] for run in runs}
        if len(samples) > 1:
            raise ValueError(f"Multiple BioSamples for {library['library_id']}")
        rows.append({
            "library_id": library["library_id"],
            "species_current": library["species_current"],
            "sra_run": library["sra_run"],
            "biosample": ";".join(sorted(samples)),
            "mapping_status": "unique_run_match" if samples else "no_public_run",
            "runinfo_sources": ";".join(sorted({source for run in runs for source in sources[run]})),
        })
    with (BASE / "biosample_mapping.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    accessions = sorted({row["biosample"] for row in rows if row["biosample"]})
    requests = []
    for start in range(0, len(accessions), 50):
        batch = accessions[start:start + 50]
        requests.append({
            "batch": f"biosamples_{start // 50 + 1:02}",
            "url": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&retmode=xml&id=" + ",".join(batch),
            "accessions": batch,
        })
    (BASE / "biosample_requests.json").write_text(json.dumps(requests, indent=2) + "\n")
    identity = {
        "prepared_utc": datetime.now(timezone.utc).isoformat(),
        "manifest_sha256": hashlib.sha256(Path("manifest.csv").read_bytes()).hexdigest(),
        "source_sha256": {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files},
    }
    (BASE / "mapping_inputs.json").write_text(json.dumps(identity, indent=2) + "\n")


if __name__ == "__main__":
    main()
