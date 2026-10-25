#!/usr/bin/env python3
"""Prepare a GTDB-to-NCBI name crosswalk using exact NCBI Taxonomy queries."""
import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


def read_table(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def request(endpoint, params):
    time.sleep(0.4)
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'+endpoint+'?'+urllib.parse.urlencode(params)
    with urllib.request.urlopen(url,timeout=60) as handle:
        return url, handle.read()


def lookup(name, cache, offline=False):
    if name in cache and (cache[name]['status'] == 'verified' or offline):
        return cache[name]
    if offline:
        value={'name':name,'ncbi_taxid':'','ncbi_scientific_name':'','rank':'',
               'query_url':'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?name='+urllib.parse.quote(name),
               'verified_utc':'','status':'pending_ncbi_verification',
               'note':'Exact-name verification pending; unavailable endpoints do not establish absence from NCBI.'}
        cache[name]=value
        return value
    url, body=request('esearch.fcgi',{'db':'taxonomy','term':name+'[Scientific Name] AND txid2[Subtree]','retmode':'json','retmax':'10'})
    payload=json.loads(body)
    ids=payload['esearchresult']['idlist']
    taxid=''
    rank=''
    scientific_name=''
    if ids:
        _, xml=request('efetch.fcgi',{'db':'taxonomy','id':','.join(ids),'retmode':'xml'})
        for taxon in ET.fromstring(xml).findall('Taxon'):
            if taxon.findtext('ScientificName','').lower()==name.lower():
                if taxid:
                    raise ValueError('Ambiguous exact NCBI name: '+name)
                taxid=taxon.findtext('TaxId','')
                rank=taxon.findtext('Rank','')
                scientific_name=taxon.findtext('ScientificName','')
    value={'name':name,'ncbi_taxid':taxid,'ncbi_scientific_name':scientific_name,'rank':rank,'query_url':url,'verified_utc':datetime.now(timezone.utc).isoformat(),'status':'verified' if taxid else 'not_found_exact_name'}
    cache[name]=value
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo',type=Path,default=Path.cwd())
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--offline', action='store_true', help='Use verified cached taxonomy records; explicitly flag all unresolved names')
    args=parser.parse_args()
    output=args.output_dir
    output.mkdir(parents=True,exist_ok=True)
    decisions=read_table(args.repo/'config/mag_deposit_decisions.tsv')
    catalog=read_table(args.repo/'data/results/post_expansion/2026-10-20/curated/phase3_catalog/mags/mag_catalog.tsv')
    catalog={r['mag_id']:r for r in catalog if r['catalog_status']=='retained'}
    names={
      'MAGSP0001':('Serratia nevei','Serratia'),
      'MAGSP0002':('Pseudomonadales bacterium','Pseudomonadales'),
      'MAGSP0003':('Planktomarina sp.','Planktomarina'),
      'MAGSP0004':('Alphaproteobacteria bacterium','Alphaproteobacteria'),
      'MAGSP0005':('Metamycoplasmataceae bacterium','Metamycoplasmataceae'),
      'MAGSP0006':('Gammaproteobacteria bacterium','Gammaproteobacteria'),
      'MAGSP0007':('Metamycoplasmataceae bacterium','Metamycoplasmataceae'),
      'MAGSP0008':('Alphaproteobacteria bacterium','Alphaproteobacteria'),
      'MAGSP0009':('Flavobacteriaceae bacterium','Flavobacteriaceae'),
      'MAGSP0010':('Metamycoplasmataceae bacterium','Metamycoplasmataceae'),
      'MAGSP0011':('Metamycoplasmataceae bacterium','Metamycoplasmataceae'),
      'MAGSP0012':('Mycoplasma sp.','Mycoplasma'),
      'MAGSP0013':('Pseudomonadales bacterium','Pseudomonadales'),
      'MAGSP0015':('Cellvibrionaceae bacterium','Cellvibrionaceae'),
      'MAGSP0016':('Crocinitomix sp.','Crocinitomix'),
      'MAGSP0017':('Cognatishimia sp.','Cognatishimia'),
      'MAGSP0019':('Pseudomonadales bacterium','Pseudomonadales'),
      'MAGSP0021':('Gammaproteobacteria bacterium','Gammaproteobacteria'),
      'MAGSP0024':('Cognatishimia sp.','Cognatishimia'),
      'MAGSP0027':('Salinisphaeraceae bacterium','Salinisphaeraceae'),
      'MAGSP0028':('Xenococcus sp.','Xenococcus'),
      'MAGSP0029':('Metamycoplasmataceae bacterium','Metamycoplasmataceae'),
      'MAGSP0030':('Pseudomonadales bacterium','Pseudomonadales'),
      'MAGSP0031':('Metamycoplasmataceae bacterium','Metamycoplasmataceae'),
      'MAGSP0033':('Hepatoplasmataceae bacterium','Hepatoplasmataceae'),
    }
    named_clusters={'MAGSP0001','MAGSP0014','MAGSP0018','MAGSP0020','MAGSP0022','MAGSP0023','MAGSP0025','MAGSP0032'}
    cache_path=output/'ncbi_taxonomy_lookup.json'
    cache=json.loads(cache_path.read_text()) if cache_path.exists() else {}
    records=[]
    for index,row in enumerate(decisions,1):
        mag=row['mag_id']; cluster=row['species_cluster']; gtdb=catalog[mag]['gtdb_classification']
        if cluster in named_clusters:
            species=next(x[3:] for x in gtdb.split(';') if x.startswith('s__'))
            proposed,parent=species,species.split()[0]
            basis='GTDB placement and reviewed type-strain ANI; PGAP taxcheck pending'
        else:
            proposed,parent=names[cluster]
            basis='Lowest proposed NCBI lineage consistent with reviewed placement; no GTDB placeholders used'
        match=lookup(proposed,cache,args.offline)
        if cluster in named_clusters and match['status']=='not_found_exact_name':
            proposed=parent+' sp.'
            match=lookup(proposed,cache,args.offline)
            basis+='; species name not found in NCBI, genus-level fallback'
        ancestor=match if match['ncbi_taxid'] else lookup(parent,cache,args.offline)
        if not ancestor['ncbi_taxid'] and not match['ncbi_taxid']:
            screening_parent = ('Alphaproteobacteria' if 'c__Alphaproteobacteria' in gtdb
                                else 'Gammaproteobacteria' if 'c__Gammaproteobacteria' in gtdb
                                else 'Mycoplasmatales' if cluster == 'MAGSP0033'
                                else 'Cyanobacteriota' if cluster == 'MAGSP0028' else '')
            if screening_parent:
                ancestor=lookup(screening_parent,cache,args.offline)
                basis+='; broader verified screening lineage used while the proposed parent is unresolved'
        status=('verified' if match['ncbi_taxid'] else
                'requires_ncbi_name_creation' if match['status']=='not_found_exact_name'
                else 'pending_ncbi_verification')
        records.append({'mag_id':mag,'isolate':f'SIPH{index:04d}','species_cluster':cluster,'gtdb_classification':gtdb,'proposed_ncbi_name':proposed,'ncbi_taxid':match['ncbi_taxid'],'name_status':status,'screening_taxid':match['ncbi_taxid'] or ancestor['ncbi_taxid'],'screening_name':match['ncbi_scientific_name'] or ancestor['ncbi_scientific_name'],'screening_taxonomy_url':match['query_url'] if match['ncbi_taxid'] else ancestor['query_url'],'basis':basis,'taxonomy_url':match['query_url'],'verified_utc':match['verified_utc']})
        cache_path.write_text(json.dumps(cache,indent=2)+'\n')
    for name, rows in [('mag_organism_names.tsv',records),('ncbi_name_verification_request.tsv',[{'isolate':r['isolate'],'proposed_organism':r['proposed_ncbi_name'],'gtdb_lineage':r['gtdb_classification'],'name_status':r['name_status']} for r in records])]:
        with (output/name).open('w') as handle:
            w=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    print(f'Wrote {len(records)} MAG crosswalk records; '+str(sum(r['name_status']=='verified' for r in records))+' have verified organism names and '+str(sum(bool(r['screening_taxid']) for r in records))+' have verified screening taxids. No email sent.')


if __name__=='__main__':
    main()
