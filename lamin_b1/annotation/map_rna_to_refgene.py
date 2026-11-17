"""Map RNA gene IDs to reciprocal-unique RefGene loci with shared exons."""
from pathlib import Path
import argparse
import collections
import csv
import sqlite3
import subprocess
import numpy as np

SAMPLES=['Basal_A','Basal_B','Basal_C','Bic_A','Bic_B','Bic_C','Bic_nimo_A','Bic_nimo_B','Bic_nimo_C']

def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle,delimiter='\t'))

def write(path, rows):
    with path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def metas(p):
 rows=read(p/'gene_metadata.tsv')
 for d in rows:
  for k in ['start0','end0']:d[k]=int(d[k])
  d['tss0']=d['start0'] if d['strand']=='+' else d['end0']-1
 assert len({d['locus_id'] for d in rows})==len(rows)
 return {d['locus_id']:d for d in rows}

def exons(p):
 out=collections.defaultdict(list)
 with (p/'exons.bed').open() as f:
  for line in f:
   x=line.rstrip().split('\t');out[x[3]].append((int(x[1]),int(x[2])))
 for k,v in out.items():
  v.sort();assert all(a<b for a,b in v) and all(v[i][1]<=v[i+1][0] for i in range(len(v)-1))
 return out

def overlap(a,b):
 i=j=total=0
 while i<len(a) and j<len(b):
  total+=max(0,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
  if a[i][1]<b[j][1]:i+=1
  else:j+=1
 return total

def joined(s):return ';'.join(sorted(map(str,s)))

def build(refgene, novogene, database, expression, output):
 ref=metas(refgene);nov=metas(novogene);re=exons(refgene);ne=exons(novogene)
 raw=read(expression);rna={d['gene_id']:d for d in raw}
 if len(rna)!=len(raw):raise ValueError('RNA gene IDs must be unique.')
 byid=collections.defaultdict(list)
 for d in nov.values():byid[d['raw_gene_id']].append(d)
 conn=sqlite3.connect('file:'+str(database)+'?mode=ro',uri=True)
 dbmeta=dict(conn.execute('select name,value from metadata'))
 def idmap(sql):
  out=collections.defaultdict(set)
  for key,entrez in conn.execute(sql):out[key].add(entrez)
  return out
 ens=idmap('select ensembl_id,gene_id from ensembl join genes using(_id)')
 acc=idmap('select accession,gene_id from refseq join genes using(_id)')
 names=idmap('select symbol,gene_id from gene_info join genes using(_id)')
 aliases=idmap('select alias_symbol,gene_id from alias join genes using(_id)')
 for k,v in aliases.items():names[k]|=v
 ref_ids={k:set().union(*(acc[t.split('.')[0]] for t in d['transcript_ids'].split(','))) for k,d in ref.items()}
 conn.close()
 # Enumerate every same-strand gene-body overlap using bedtools.
 command=['bedtools','intersect','-a',str(novogene/'genes.bed'),'-b',str(refgene/'genes.bed'),'-s','-wo']
 proc=subprocess.run(command,text=True,capture_output=True,check=True)
 pairs={}
 for line in proc.stdout.splitlines():
  x=line.split('\t');assert len(x)==15;p=(x[3],x[10]);assert p not in pairs;pairs[p]=int(x[14])
 evidence=[];eligible=collections.defaultdict(set);reverse=collections.defaultdict(set);pair_by={}
 for (nid,rid),bp in sorted(pairs.items()):
  n=nov[nid];r=ref[rid];gid=n['raw_gene_id'];same=n['gene_name']==r['gene_name']
  eni=ens[gid];rsi=ref_ids[rid];ei=eni & rsi
  alias=(eni or names[n['gene_name']]) & (rsi or names[r['gene_name']])
  conflict=bool(eni and rsi and not ei)
  shared=overlap(ne[nid],re[rid]);identity=same or bool(ei) or bool(alias)
  qualifies=gid in rna and len(byid[gid])==1 and identity and shared>0 and not conflict
  method='Ensembl_RefSeq_Entrez' if ei else ('exact_symbol' if same else ('local_alias_Entrez' if alias else 'none'))
  d=dict(rna_gene_id=gid,novogene_locus_id=nid,rna_name=n['gene_name'],refgene_locus_id=rid,refgene_name=r['gene_name'],chromosome=n['chrom'],strand=n['strand'],gene_overlap_bp=bp,shared_exon_bp=shared,exon_fraction_novogene=shared/sum(b-a for a,b in ne[nid]),exon_fraction_refgene=shared/sum(b-a for a,b in re[rid]),exact_symbol=same,stable_ID_conflict=conflict,rna_Entrez=joined(eni),refseq_Entrez=joined(rsi),common_Entrez=joined(ei or alias),identity_evidence=method,identity_supported=identity,eligible_before_uniqueness=qualifies)
  evidence.append(d);pair_by[nid,rid]=d
  if qualifies:eligible[gid].add(rid);reverse[rid].add(gid)
 matched={gid:next(iter(rs)) for gid,rs in eligible.items() if len(rs)==1 and len(reverse[next(iter(rs))])==1}
 assert len(matched)==len(set(matched.values()))
 statuses={};cross=[];impact=[]
 bypair=collections.defaultdict(list)
 for d in evidence:bypair[d['rna_gene_id']].append(d)
 for gid,d in rna.items():
  if gid in matched:status='matched_one_to_one'
  elif not byid[gid]:status='no_novogene_locus'
  elif len(byid[gid])>1:status='multiple_novogene_loci'
  elif len(eligible[gid])>1:status='multiple_refgene_candidates'
  elif len(eligible[gid])==1:status='shared_refgene_candidate'
  elif any(e['identity_supported'] and e['stable_ID_conflict'] for e in bypair[gid]):status='identity_conflict'
  elif any(e['identity_supported'] and e['shared_exon_bp']==0 for e in bypair[gid]):status='no_shared_exons'
  else:status='no_verified_refgene_match'
  statuses[gid]=status
  row=dict(rna_gene_id=gid,rna_name=d['gene_name'],gene_biotype=d['gene_biotype'],status=status,novogene_loci=joined(x['locus_id'] for x in byid[gid]),eligible_refgene_candidates=joined(eligible[gid]),refgene_locus_id=matched.get(gid,''),refgene_name='',chromosome='',strand='',refgene_start0='',refgene_end0='',refgene_tss0='',evidence='',shared_exon_bp='',exon_fraction_novogene='',exon_fraction_refgene='')
  if gid in matched:
   rid=matched[gid];r=ref[rid];n=byid[gid][0];e=pair_by[n['locus_id'],rid]
   row.update(refgene_name=r['gene_name'],chromosome=r['chrom'],strand=r['strand'],refgene_start0=r['start0'],refgene_end0=r['end0'],refgene_tss0=r['tss0'],evidence=e['identity_evidence'],shared_exon_bp=e['shared_exon_bp'],exon_fraction_novogene=e['exon_fraction_novogene'],exon_fraction_refgene=e['exon_fraction_refgene'])
  cross.append(row)
 output.mkdir(parents=True,exist_ok=True)
 write(output/'rna_refgene_crosswalk.tsv',cross)
 genes=[]
 for row in cross:
  raw=rna[row['rna_gene_id']]
  values=np.array([float(raw[name]) for name in SAMPLES])
  if not (np.isfinite(values).all() and (values>=0).all()):
   raise ValueError('FPKM values must be finite and nonnegative.')
  d=dict(gene_id=row['rna_gene_id'],gene_name=row['rna_name'],gene_biotype=row['gene_biotype'],
         annotation_status='unique' if row['status']=='matched_one_to_one' else row['status'],
         refgene_locus_id=row['refgene_locus_id'],refgene_name=row['refgene_name'],mapping_evidence=row['evidence'],
         chromosome=row['chromosome'],start0=row['refgene_start0'],end0=row['refgene_end0'],strand=row['strand'],tss0=row['refgene_tss0'])
  d.update({name+'_FPKM':float(value) for name,value in zip(SAMPLES,values)})
  for i,name in enumerate(['TTX','Bic','Bic_Nimo']):d[name+'_mean_FPKM']=float(values[i*3:i*3+3].mean())
  genes.append(d)
 write(output/'gene_coordinates_expression.tsv',genes)
 print(f'Mapped {len(matched):,} RNA genes to reciprocal-unique RefGene loci.')

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 for name in ['refgene','novogene','database','expression','output']:
  parser.add_argument('--'+name,type=Path,required=True)
 args=parser.parse_args()
 build(args.refgene,args.novogene,args.database,args.expression,args.output)

if __name__=='__main__':
 main()
