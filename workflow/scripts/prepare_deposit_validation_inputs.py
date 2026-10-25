#!/usr/bin/env python3
"""Prepare draft table2asn inputs without inventing submitter or feature metadata."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from package_deposit_sequences import write_fasta, write_table
from prepare_deposit_mag import read_fasta, read_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle
    rows = read_table(bundle / 'mag_manifest.tsv')
    metadata = {r['isolate']: r for r in read_table(bundle / 'mimag_biosample_metadata.tsv')}
    destination = bundle / 'validation_inputs'
    manifest = []
    for row in rows:
        if row['route'] == 'held_provenance_only':
            manifest.append({'unit': row['mag_id'], 'route': row['route'], 'fasta': row['fasta'],
                             'source_table': '', 'feature_table': '',
                             'status': 'held_provenance_only_not_submission_or_claim_archive',
                             'source_fields_missing': 'sequence_identity_or_structure_requires_review',
                             'pending': row['blockers']})
            continue
        source = metadata[row['isolate']]
        sequences = read_fasta(Path(row['fasta']))
        base = destination / row['route'] / row['isolate']
        fasta = base.with_suffix('.fsa')
        write_fasta(fasta, list(sequences.items()))
        sources = [{'SeqID': name, 'organism': source['organism'], 'isolate': source['isolate'],
                    'host': source['host'], 'collection_date': source['collection_date'],
                    'geo_loc_name': source['geo_loc_name'], 'lat_lon': source['lat_lon'],
                    'isolation_source': source['isolation_source']} for name in sequences]
        source_table = base.with_suffix('.src')
        write_table(source_table, sources)
        missing = [field for field in ('organism', 'collection_date', 'geo_loc_name', 'isolation_source')
                   if not source[field]]
        manifest.append({'unit': row['mag_id'], 'route': row['route'], 'fasta': str(fasta),
                         'source_table': str(source_table), 'feature_table': '',
                         'status': 'draft_not_validated', 'source_fields_missing': ';'.join(missing),
                         'pending': row['blockers'] + ';real_submitter_template_sbt_required'})
    for row in read_table(bundle / 'sequence_sets.tsv'):
        if row['set'] == 'mags':
            continue
        fasta = destination / 'other_sequences' / f"{row['set']}.fsa"
        write_fasta(fasta, list(read_fasta(Path(row['fasta'])).items()))
        manifest.append({'unit': row['set'], 'route': 'sequence_specific_review', 'fasta': str(fasta),
                         'source_table': '', 'feature_table': '', 'status': 'draft_not_validated',
                         'source_fields_missing': 'NCBI_organism_and_per_sequence_source_qualifiers',
                         'pending': 'source_evidence_review;feature_coordinates_and_translation_review;'
                                    'real_submitter_template_sbt_required;table2asn;'
                                    + ('ribovore' if row['set'] == 'ssu' else 'sequence_route_review')})
    write_table(destination / 'validation_inputs.tsv', manifest)
    status = {
        'status': 'draft_inputs_only_not_validated',
        'table2asn_template': 'pending_real_author_submitter_template_no_template_invented',
        'external_validation_tools': ['table2asn', 'PGAP', 'ribovore'],
        'tool_execution_status': 'not_assessed_by_draft_input_preparation',
        'MAG_draft_inputs': sum(r['route'] != 'held_provenance_only' for r in rows),
        'held_provenance_only': sum(r['route'] == 'held_provenance_only' for r in rows),
        'non_MAG_feature_tables': 'pending_review_of_preserved_source_evidence',
        'next': 'Review FCS reports, install validation tools, and supply real template before table2asn.',
        'format_reference': 'https://www.ncbi.nlm.nih.gov/genbank/table2asn/',
    }
    (destination / 'status.json').write_text(json.dumps(status, indent=2) + '\n')


if __name__ == '__main__':
    main()
