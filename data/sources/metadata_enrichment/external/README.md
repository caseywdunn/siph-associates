# Public collection metadata retrieval

BioSample records were retrieved on 2026-10-04 to supplement the specimen
spreadsheets. These files supply candidates for review; they do not change
the project manifest or any biological result on their own. The GBIF
subdirectory is a separate retrieval maintained alongside these records.

## Inputs and retrieval

`biosample_mapping.tsv` joins the 205 current library IDs to exact public run
accessions in `manifest.csv` and the retained, headerless RunInfo snapshots
in `data/sources/*runinfo*`. It identifies 185 unique BioSamples for 185
libraries. Twenty libraries lack a public run accession and were not given
an inferred BioSample ID. The 21 libraries noted as unreleased in the
historical Church study differ from these 20 gaps because YPM-IZ-104465 has
a reused public Ahuja run, SRR23143271. The empty Church full-RunInfo file
supplies no rows; other retained snapshots supply the validated mappings.

`mapping_inputs.json` records the input hashes. The mapping and exact NCBI
E-utilities URLs can be recreated from the analysis repository root:

```bash
python3 data/sources/metadata_enrichment/external/prepare_biosample_requests.py
```

The four URLs in `biosample_requests.json` were fetched with
`curl --fail --location --connect-timeout 15 --max-time 45 URL -o FILE`,
writing `biosamples_01.xml` through `biosamples_04.xml`. All four requests
completed successfully. The raw priority check is also retained at
`raw/biosamples_priority.xml`; its four accessions are included in the
complete batches, and it is not counted as additional records.

Normalize the complete retained XML snapshots offline:

```bash
python3 data/sources/metadata_enrichment/external/normalize_biosamples.py
```

Both scripts use the Python standard library. The normalizer validates each
batch's full returned accession set against its request before emitting any
candidates. `biosample_provenance.json` records exact URLs, raw checksums,
sizes, retrieval times from file modification timestamps, code identity and
coverage. It also records the counts of the emitted fields.

## Interpretation and handoff

The review input is `biosample_candidates.tsv`; raw attribute names and
values are retained in `biosample_attributes.tsv`, and retrieval coverage is
in `biosample_coverage.tsv`. `biosample_issues.tsv` records eight specific
issues without altering the reported values. Source URLs and attribute
locators accompany each candidate. No host taxonomy is reassigned from
BioSample organism names.

There are 35 explicit metre-valued depth points and 14 explicit metre
intervals. Six records provide bare depth numbers without units: WS5,
WS6, WS7, WS8, WS9 and WS10. These supply `depth_original` and provenance
only, with no numeric depth inferred. In particular, the unusually deep
WS5 and WS6 values repeat their spreadsheet submissions and do not resolve
whether those numbers denote collection or seafloor depth. NA19's BioSample
reports 297 metres, agreeing with the 2024 sheet but conflicting with the
2026 value of 296. CWD16's 10-metre point lies within the 2026 0–20-metre
range; that is a compatible pair of descriptions, not a point disagreement.

The Guam specimen YPM-IZ-104465 repeats the submitted west longitude
`13.428 N 144.799 W`. It is retained and flagged, never silently corrected
to east. Missing-value strings are excluded from candidates. Coordinates
are converted only by applying the explicit hemisphere signs. Original
`geo_loc_name` strings populate locality, with submitted country/territory
prefixes also retained; they do not alter the study's `ocean_region`
inference categories. `isolation_source` is retained as a raw attribute,
not interpreted as a sampling method.

The records provide 185 collection dates, 179 coordinate pairs and 169
nonmissing tissues. No beach collection method is documented by these
BioSample attributes: `isolation_source` is either `ocean` or `not
applicable`. Therefore no new blanket zero-depth assumption is made for
Physalia. Metadata conflict review and promotion are handled by the parent
metadata-enrichment workflow.
