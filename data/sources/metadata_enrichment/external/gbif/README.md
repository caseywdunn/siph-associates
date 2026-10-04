# GBIF exact-voucher metadata check

On 2026-10-04, the GBIF occurrence API was searched for the cohort's exact
Yale Peabody Museum invertebrate zoology vouchers in dataset
`854e35e6-f762-11e1-a439-00145eb45e9a`, **Invertebrate Zoology Division,
Yale Peabody Museum** ([dataset DOI](https://doi.org/10.15468/0lkr3w)).
This was a bounded collection-metadata check; no sequence analyses were run.

## Coverage and result

The cohort contains 205 libraries. Six successful batches searched the
168 libraries with exact YPM IZ voucher identifiers. Each batch returned
`count: 0`, an empty result list and `endOfRecords: true`. Thus no cohort
voucher matched the published dataset. The remaining 37 libraries lack an
exact YPM IZ voucher and were not searched in this dataset. They are
explicitly distinguished from successful searches with no match in
`gbif_coverage.tsv`; there were no failed cohort requests.

Catalog strings were queried in the observed GBIF format, `YPM IZ 111736`.
Both five-digit and zero-padded six-digit forms were included when relevant.
A positive control, `YPM IZ 116025`, returned one occurrence. A second
control with both `YPM IZ 116025` and `YPM IZ 116024` returned two, confirming
that repeated `catalogNumber` parameters use OR matching. The dataset
sample returned three records from a total of 150,964 published occurrences.

A supplemental dataset-restricted `q=Physalia` search returned 80 records
and reached the end of the result set; none of their catalog numbers
matched the cohort. The equivalent Nanomia search returned zero records.
Neither taxon-only search was used to assign specimen metadata.

`gbif_candidates.tsv` deliberately contains its header only: GBIF supplies
no verified metadata additions for this cohort. These negative searches do
not establish that the specimens or museum collection records do not exist;
they describe the GBIF-published dataset at retrieval. In particular, no
beach/surface collection designation or zero-metre depth can be inferred
from unrelated Physalia occurrences.

## Provenance

`cohort_requests.json` preserves the six exact cohort URLs and the first
control URL, with their canonical library IDs. Raw JSON responses are kept
alongside this file. `gbif_provenance.json` records every retained query,
response checksum, byte count, retrieval date, filesystem timestamp, result
count and control outcome. Filesystem timestamps are not HTTP response
timestamps; HTTP response headers were not retained. Supplemental request
URLs are recorded from the retrieval task, with their provenance identified.

`gbif_coverage.tsv` has one row per canonical library, including the searched
catalog variants, exact batch URL, response file and checksum. Its manifest
checksum is recorded in `gbif_provenance.json`. Catalog comparison permits
only the recognized YPM/IZ separators and leading zeros; host names and
partial identifiers are never used to assign metadata. The shared candidate
schema is `library_id`, `field`, `value`, `source_type`, `source`,
`source_locator`, `source_record`, `notes`.
