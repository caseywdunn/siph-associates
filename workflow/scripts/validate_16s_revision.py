#!/usr/bin/env python3
"""Verify the sole reviewed 16S query deletion and the reported sample totals."""
import argparse
import csv
import json
from pathlib import Path

from prepare_deposit_mag import read_fasta


def sample_counts(path):
    with path.open(newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    return {'samples':len({r['sample'] for r in rows}),
            'publications':len({r['Publication'] for r in rows}),
            'matching_rows':len(rows)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--before',required=True,type=Path)
    parser.add_argument('--after',required=True,type=Path)
    parser.add_argument('--report',required=True,type=Path)
    args=parser.parse_args()
    old=read_fasta(args.before/'mycoplasmatales_mag_16s.fasta')
    new=read_fasta(args.after/'mycoplasmatales_mag_16s.fasta')
    expected={name for name in old if '|Church2025__WAM-Z97916__bin.2|' in name
              and name.split('|')[-1].split(':')[0].split('__')[-1]=='k141_5955'}
    if len(expected)!=1 or set(old)-set(new)!=expected or set(new)-set(old):
        raise ValueError(f'Unexpected 16S query change: expected={expected}; removed={set(old)-set(new)}; added={set(new)-set(old)}')
    if any(sequence != old[name] for name,sequence in new.items()):
        raise ValueError('A retained 16S sequence changed during the query correction')
    before=sample_counts(args.before/'matching_asv_samples.tsv')
    after=sample_counts(args.after/'matching_asv_samples.tsv')
    if any(counts['samples']!=98 or counts['publications']!=23 for counts in [before,after]):
        raise ValueError(f'Expected 98 samples and 23 publications in both runs: before={before}; after={after}')
    def attributions(folder):
        with (folder/'matching_asv_samples.tsv').open(newline='') as handle:
            rows=list(csv.DictReader(handle,delimiter='\t'))
        result={(r['library'],r['sample'],r['asv']):r for r in rows}
        if len(result)!=len(rows):
            raise ValueError('Duplicate library/sample/ASV rows in matching table')
        return result
    old_rows,new_rows=attributions(args.before),attributions(args.after)
    changed=[]
    for key in sorted(old_rows.keys() & new_rows.keys()):
        a,b=old_rows[key],new_rows[key]
        if a['species_cluster']!=b['species_cluster']:
            changed.append({'library':key[0],'sample':key[1],'asv':key[2],
                            'host':b['Scientific_Name'],'before':a['species_cluster'],
                            'after':b['species_cluster']})
    report={'removed_query':sorted(expected),'queries_before':len(old),'queries_after':len(new),
            'before':before,'after':after,'retained_sequences_unchanged':True,
            'changed_attributions':changed,
            'removed_matching_rows':len(old_rows.keys()-new_rows.keys()),
            'added_matching_rows':len(new_rows.keys()-old_rows.keys()),
            'verified_scope':'misbinned clade-B query deletion; headline sample/publication counts'}
    args.report.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
