#!/usr/bin/env python3
"""Export the corrected MAGSP0031 size and freshly measured CheckM2 quality."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

from prepare_deposit_mag import read_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path('data/results/deposit'))
    parser.add_argument('--output', type=Path, default=Path('tables/deposit/magsp0031_revision.tsv'))
    parser.add_argument('--execution', required=True)
    args = parser.parse_args()
    mag = 'Church2025__WAM-Z97910__bin.1'
    edit_path = args.root / 'edits' / f'{mag}.json'
    edit = json.loads(edit_path.read_text())
    quality_path = args.root / 'checkm2' / mag / 'quality_report.tsv'
    quality = read_table(quality_path)
    if len(quality) != 1 or quality[0]['Name'] != mag:
        raise ValueError('CheckM2 report does not identify the corrected MAGSP0031')
    if int(quality[0]['Genome_Size']) != edit['genome_size']:
        raise ValueError('CheckM2 genome size differs from the corrected FASTA')
    version = (Path(args.execution) / 'checkm2_version.txt').read_text().splitlines()[-1]
    if not re.fullmatch(r'\d+\.\d+\.\d+', version):
        raise ValueError('Unrecognized recorded CheckM2 version')
    fasta = args.root / 'cleaned_mags' / f'{mag}.fa'
    if hashlib.sha256(fasta.read_bytes()).hexdigest() != edit['output_sha256']:
        raise ValueError('Corrected FASTA does not match the sequence-edit record')
    row = {
        'mag_id': mag, 'species_cluster': 'MAGSP0031',
        'original_genome_size': edit['original_size'],
        'corrected_genome_size': edit['genome_size'],
        'removed_duplicate_bases': edit['original_size'] - edit['genome_size'],
        'contigs': edit['contigs'], 'completeness': quality[0]['Completeness'],
        'contamination': quality[0]['Contamination'],
        'checkm2_version': version, 'checkm2_report': str(quality_path),
        'checkm2_report_sha256': hashlib.sha256(quality_path.read_bytes()).hexdigest(),
        'sequence_edit_record': str(edit_path), 'corrected_fasta': str(fasta),
        'original_fasta_sha256': edit['input_sha256'],
        'corrected_fasta_sha256': edit['output_sha256'],
        'execution': args.execution, 'mapping_catalog_unchanged': 'true',
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerow(row)


if __name__ == '__main__':
    main()
