"""Extract public cortical CpG counts with the published coordinate and coverage rules."""
from pathlib import Path
import argparse,gzip,collections,json
import numpy as np,pandas as pd
import subprocess,pysam
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs',type=Path,required=True)
parser.add_argument('--public',type=Path,required=True)
parser.add_argument('--fasta',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--reference-mask',type=Path)
parser.add_argument('--chain',type=Path)
args=parser.parse_args();inputs=args.inputs;P=args.public;R=args.output;R.mkdir(parents=True,exist_ok=True)
chroms=['chr'+str(i) for i in range(1,20)]
mask=np.load(args.reference_mask or inputs/'reference_mask.npz');sizes=dict(zip(mask['chromosomes'],mask['lengths']));masks={}
for c in chroms:
 a=[]
 for st,en in sorted((int(st),int(en)) for ty in ['mask','nonACGT'] for st,en in zip(mask[c+'_'+ty+'_start'],mask[c+'_'+ty+'_end'])):
  if a and st<=a[-1][1]:a[-1][1]=max(en,a[-1][1])
  else:a.append([st,en])
 masks[c]=np.asarray(a,dtype=np.int64)
blocks=collections.defaultdict(list);h=None
for line in gzip.open(args.chain or inputs/'mm9ToMm10.over.chain.gz','rt'):
 a=line.split()
 if not a:continue
 if a[0]=='chain':h=a;t=int(a[5]);q=int(a[10]);continue
 z=int(a[0])
 if h[2] in chroms and h[7] in chroms:
  assert h[4]=='+'
  blocks[h[2]].append((t,t+z,chroms.index(h[7]),q,int(h[8]),h[9]=='-'))
 t+=z;q+=z
 if len(a)==3:t+=int(a[1]);q+=int(a[2])
for c in blocks:
 blocks[c]=np.asarray(sorted(blocks[c]),dtype=np.int64);assert np.all(blocks[c][1:,0]>=blocks[c][:-1,1])
fa=pysam.FastaFile(str(args.fasta));seq_cache={}
def seq(c):
 if c not in seq_cache:
  seq_cache.clear();seq_cache[c]=np.frombuffer(fa.fetch(c).upper().encode(),dtype=np.uint8)
 return seq_cache[c]
def is_cpg(c,pos):
 s=seq(c);valid=(pos>=1)&(pos<len(s)-1);p=np.clip(pos,1,len(s)-2)
 return valid&(((s[p]==67)&(s[p+1]==71))|((s[p]==71)&(s[p-1]==67)))
for age,sample in [('1wk','GSM2757470'),('2wk','GSM2757471')]:
 report=R/f'public_{age}_summary.json'
 if report.exists():continue
 source=P/f'{sample}_Cortex_BS_{age}.txt.gz';assert source.exists(),source
 cache=R/f'public_{age}_CpG_mm10';cache.mkdir(exist_ok=True)
 qc=collections.Counter();hist=np.zeros(100000,dtype=np.int64);parts=collections.defaultdict(list)
 if not (cache/'extraction.json').exists():
  print('Extracting public',age,'CG records and converting coordinates',flush=True)
  err=(R/f'{age}_extract.log').open('wb')
  p1=subprocess.Popen(['gzip','-dc',str(source)],stdout=subprocess.PIPE,stderr=err)
  p2=subprocess.Popen(['awk','$3=="CG"'],stdin=p1.stdout,stdout=subprocess.PIPE,stderr=err);p1.stdout.close()
  for d in pd.read_csv(p2.stdout,sep='\t',header=None,names=['chr','pos','context','meth','depth'],dtype={'chr':'category','pos':'int64','context':'category','meth':'int64','depth':'int64'},chunksize=500000):
   qc['CG_records']+=len(d);assert (d.meth<=d.depth).all();assert (d.meth>=0).all()
   for c,g in d.groupby('chr',observed=True):
    if c not in blocks:continue
    qc['autosomal_CG']+=len(g);pos=g.pos.to_numpy()-1;b=blocks[c];ix=np.searchsorted(b[:,0],pos,side='right')-1;k=np.maximum(ix,0);good=(ix>=0)&(pos<b[k,1]);qc['source_unmapped']+=int((~good).sum())
    a=b[k[good]];dist=pos[good]-a[:,0];p=np.where(a[:,5]>0,a[:,4]-1-(a[:,3]+dist),a[:,3]+dist);meth=g.meth.to_numpy()[good];depth=g.depth.to_numpy()[good]
    for code in np.unique(a[:,2]):
     take=a[:,2]==code;target=chroms[int(code)];pp=p[take];mm=meth[take];nn=depth[take]
     valid=is_cpg(target,pp);qc['lifted_records']+=len(pp);qc['mm10_CpG_valid']+=int(valid.sum());qc['mm10_context_changed']+=int((~valid).sum())
     pp,mm,nn=pp[valid],mm[valid],nn[valid]
     keep=(pp//5000<int(sizes[target])//5000);m=masks[target];mi=np.searchsorted(m[:,0],pp,side='right')-1;keep&=(mi<0)|(pp>=m[np.maximum(mi,0),1]);qc['masked_or_partial_window']+=int((~keep).sum())
     pp,mm,nn=pp[keep],mm[keep],nn[keep];assert nn.max(initial=0)<len(hist);hist+=np.bincount(nn,minlength=len(hist))
     parts[target].append(np.column_stack((pp,mm,nn)).astype(np.int32))
   if qc['CG_records']%5000000==0:print(age,qc['CG_records'],'CpG records processed',flush=True)
  p2.stdout.close();codes=[p2.wait(),p1.wait()];err.close();assert codes==[0,0],codes
  assert qc['mm10_CpG_valid']/qc['lifted_records']>.98,dict(qc)
  for c in chroms:
   arr=np.concatenate(parts.pop(c));arr=arr[np.argsort(arr[:,0],kind='stable')]
   _,inv,multiplicity=np.unique(arr[:,0],return_inverse=True,return_counts=True);unique=multiplicity[inv]==1
   qc['ambiguous_destination_records_excluded']+=int((~unique).sum());arr=arr[unique]
   assert np.all(np.diff(arr[:,0])>0)
   np.savez_compressed(cache/f'{c}.npz',pos0=arr[:,0],meth=arr[:,1],depth=arr[:,2])
  np.save(cache/'depth_hist.npy',hist);(cache/'extraction.json').write_text(json.dumps(dict(qc),indent=2)+'\n')
 else:qc=collections.Counter(json.loads((cache/'extraction.json').read_text()));hist=np.load(cache/'depth_hist.npy')
 hist[:4]=0;cutoff=int(np.searchsorted(np.cumsum(hist),hist.sum()*.999));qc['upper_coverage_exclusive']=cutoff
 report.write_text(json.dumps(dict(qc),indent=2)+'\n')
