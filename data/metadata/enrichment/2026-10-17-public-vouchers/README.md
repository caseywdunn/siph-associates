# Public voucher recheck, October 17, 2026

Fresh requests to public LUX and GBIF APIs again found 148 of the 168 YPM vouchers
in LUX and none in GBIF. The same 20 vouchers listed in the October 12 audit
remain unconfirmed. The requested public release is not yet visible in these
search results; this does not establish the status of museum internal records.

`voucher_public_access.tsv` has exact voucher matches and public URLs;
`unconfirmed_vouchers.tsv` has the unresolved list. Both services' known-public
positive controls passed. Searches followed all pages and inspected candidate
records for exact catalog identifiers. The LUX fallback also searched outside
the Invertebrate Zoology collection restriction. Request URLs, retrieval times,
response checksums, and raw JSON are retained in `requests.json` and `raw/`.

The audit ran before applying the October 17 metadata additions. Its input
checksum matches `../2026-10-17-lux/sample_metadata_before.tsv`; voucher and
library IDs were unchanged by those additions.

```bash
python3 scripts/check_public_ypm_vouchers.py \
  --output-dir data/metadata/enrichment/2026-10-17-public-vouchers
```

Use a new dated directory for another fresh check. No earlier responses were
reused for this recheck.
