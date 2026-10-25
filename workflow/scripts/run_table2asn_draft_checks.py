#!/usr/bin/env python3
"""Run table2asn on candidate genomes, retaining missing-template status explicitly."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

from package_deposit_sequences import write_table
from prepare_deposit_mag import read_table


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, required=True)
    parser.add_argument('--table2asn', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--submitter-template', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.submitter_template and not args.submitter_template.is_file():
        raise ValueError('The supplied author submission template does not exist')
    rows = []
    for candidate in read_table(args.bundle / 'genbank_candidates.tsv'):
        isolate = candidate['isolate']
        source = args.bundle / 'validation_inputs/genbank_candidate' / isolate
        destination = args.output / isolate
        destination.mkdir()
        for suffix in ('.fsa', '.src'):
            shutil.copy2(source.with_suffix(suffix), destination / f'{isolate}{suffix}')
        fasta = destination / f'{isolate}.fsa'
        command = [str(args.table2asn), '-i', str(fasta), '-outdir', str(destination),
                   '-M', 'n', '-V', 'vb', '-Z']
        if args.submitter_template:
            command += ['-t', str(args.submitter_template)]
        (destination / 'command.json').write_text(json.dumps(command, indent=2) + '\n')
        with (destination / 'table2asn.log').open('w') as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
        reports = sorted(destination.glob('*.val'))
        stats = sorted(destination.glob('*.stats'))
        if result.returncode != 0 or not (reports or stats):
            (args.output / 'runtime_failure.json').write_text(json.dumps({
                'mag_id': candidate['mag_id'], 'exit_code': result.returncode,
                'submitter_template_supplied': bool(args.submitter_template),
                'log': str(destination / 'table2asn.log'),
                'status': 'tool_execution_incomplete_not_validated'}, indent=2) + '\n')
            raise RuntimeError(f'table2asn did not complete a validation report for {isolate}; inspect its log')
        text = '\n'.join(p.read_text() for p in reports)
        rows.append({'mag_id': candidate['mag_id'], 'isolate': isolate,
                     'input_fasta_sha256': hashlib.sha256(fasta.read_bytes()).hexdigest(),
                     'tool_exit_code': result.returncode,
                     'error_or_reject_lines': sum(bool(re.search(r'\b(ERROR|REJECT|FATAL)\b', line, re.I))
                                                  for line in text.splitlines()),
                     'warning_lines': sum('WARNING' in line.upper() for line in text.splitlines()),
                     'reports': ';'.join(map(str, reports + stats)),
                     'status': 'reports_require_review' if args.submitter_template else 'structural_draft_no_submitter_template',
                     'submission_ready': 'false'})
        write_table(args.output / 'validation_summary.tsv', rows)
    (args.output / 'status.json').write_text(json.dumps({
        'genomes_checked': len(rows), 'submission_ready': False,
        'submitter_template_supplied': bool(args.submitter_template),
        'interpretation': 'Tool reports are retained for structural and source-field review; no submission-ready status is inferred.',
        'non_MAG_source_and_feature_validation': 'pending_sequence_specific_annotations',
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
