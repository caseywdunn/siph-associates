# Public YPM voucher records

This audit checks the 168 distinct YPM Invertebrate Zoology vouchers in
`figures/cohort/sample_metadata.tsv` against public Yale LUX and GBIF records.
It reads the metadata table and does not run sequence analyses or change their
inputs. Retrieval dates, request URLs, response checksums, and input and script
checksums are retained with the results.

Read `voucher_public_access.tsv` for the result and public record links for each
voucher. `unconfirmed_vouchers.tsv` lists only vouchers without an exact match
in either service. Counts and interpretation are in `summary.json`.

A match requires the exact institution, collection, and catalog number in the
returned specimen record. Spaces, colons, hyphens, and leading zeroes are
normalized; taxon names and numeric search hits alone do not establish identity.
Museum taxon labels are preserved separately from the study's current taxonomy.

LUX searches use both padded and unpadded catalog numbers, initially restricted
to the Peabody Invertebrate Zoology collection. Unmatched vouchers are then
searched across all LUX collections. All result pages are followed, and each
candidate's full record is inspected for a matching `Identifier`. GBIF searches
use institution `YPM`, collection `IZ`, and catalog-number variants across all
datasets. Known public specimens are included as positive controls for both
services. A negative search does not establish that a specimen record is private
or absent from the museum's internal catalog.

The LUX API usage follows the public implementations in
[luxy](https://github.com/project-lux/luxy/tree/d18146f2d250a05a691039e5f84af8397fbf4af3)
and [luxmcp](https://github.com/caseywdunn/luxmcp/tree/07a4001a6c43363315906011be6380ec676d8428).
The script calls the public API directly and requires no credentials or additional
Python packages. It was run with Python 3.11.13:

```bash
python3 scripts/check_public_ypm_vouchers.py \
  --output-dir data/metadata/enrichment/2026-10-12-public-vouchers
```

Use a new dated output directory for a fresh audit. `--resume` reuses previously
retrieved responses only when their stored checksums match; it was used to add
the broader LUX search to this audit. `requests.json` records the actual retrieval
time of each response, including reused responses. Raw JSON responses are in
`raw/` under filenames derived from the request URL's SHA-256 digest.
