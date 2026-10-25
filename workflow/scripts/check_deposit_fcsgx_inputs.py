#!/usr/bin/env python3
"""Check final MAG hashes, adaptor reports and FCS-GX taxonomy divisions."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path

from package_deposit_sequences import write_table
from prepare_deposit_mag import read_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    root = Path(config['output_root'])
    names = {r['mag_id']: r for r in read_table(Path(config['organism_names']))}
    overrides = {r['ncbi_taxid']: r for r in read_table(Path(config['fcs_gx_division_overrides']))}
    manifest = read_table(root / 'bundle/mag_manifest.tsv')
    if {r['mag_id'] for r in manifest} != set(names):
        raise ValueError('Packaged MAGs and verified taxonomy crosswalk differ')
    adaptor = root / 'fcs_adaptor/mags/fcs_adaptor_report.txt'
    with adaptor.open() as handle:
        if any(line.strip() and not line.startswith('#') for line in handle):
            raise ValueError('MAG adaptor findings require review before FCS-GX')
    taxids = {r['screening_taxid'] for r in names.values()}
    divisions = {}
    with gzip.open(Path(config['databases']['fcs_gx']) / 'all.blast_div.tsv.gz', 'rt') as handle:
        for line in handle:
            fields = line.rstrip().split('\t')
            if fields[0] in taxids:
                divisions[fields[0]] = fields[1]
    rows = []
    for mag in manifest:
        name = names[mag['mag_id']]
        taxid = name['screening_taxid']
        if not taxid.isdigit() or taxid not in divisions:
            raise ValueError(f"No verified, locally supported screening taxid: {mag['mag_id']}")
        override = overrides.get(taxid)
        if divisions[taxid] == 'bacteria' and not override:
            raise ValueError(f'Unresolved broad bacterial screening division: {taxid}')
        if override and divisions[taxid] != override['database_blast_division']:
            raise ValueError(f'GX division metadata changed; review override for {taxid}')
        fasta = Path(mag['fasta'])
        if hashlib.sha256(fasta.read_bytes()).hexdigest() != mag['sha256']:
            raise ValueError(f'Packaged MAG checksum changed: {fasta}')
        rows.append({'mag_id': mag['mag_id'], 'route': mag['route'],
                     'screening_taxid': taxid, 'screening_name': name['screening_name'],
                     'database_blast_division': divisions[taxid],
                     'explicit_gx_division': override['asserted_gx_division'] if override else '',
                     'fasta_sha256': mag['sha256'], 'adaptor_review': 'no_reported_findings',
                     'status': 'inputs_ready_for_fcs_gx_not_yet_screened'})
    write_table(args.output, rows)
    print(f'FCS-GX input checks passed for {len(rows)} MAGs, including archive genomes.')


if __name__ == '__main__':
    main()
