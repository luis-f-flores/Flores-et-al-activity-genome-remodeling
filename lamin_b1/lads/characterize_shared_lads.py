"""Summarize sequence composition and gene density in shared and variable LADs."""
from pathlib import Path
import argparse
import csv
import gzip
import numpy as np
import pandas as pd
G=[[0,1],[2,3,4],[5,6,7]]

def context(signals, decoded, fasta, maskp, rmsk, tssp, output):
 import pysam
 a=dict(np.load(signals));st=np.load(decoded)['states'];sig=(np.column_stack([(st[:,g]>=1).any(1) for g in G])*[1,2,4]).sum(1);cls=np.where(sig==7,2,np.where(sig==0,0,1));cut=np.r_[True,(a['chromosome'][1:]!=a['chromosome'][:-1])|(np.diff(a['start'])!=30000)|(np.diff(cls)!=0)];owner=np.cumsum(cut)-1;beg=np.flatnonzero(cut);end=np.r_[beg[1:],len(owner)];seg=pd.DataFrame(dict(chromosome=a['chromosome'][beg],start0=a['start'][beg],end0=a['start'][end-1]+30000,class_code=cls[beg],bins=end-beg))
 seg['class']=seg.class_code.map({0:'Outside all LADs',1:'Variable LADs',2:'Shared LADs'})
 masks={c:[] for c in np.unique(a['chromosome'])};lines={c:[] for c in masks};tss={c:[] for c in masks}
 for r in csv.reader(maskp.open(),delimiter='\t'):
  if r[0] in masks:masks[r[0]].append((int(r[1]),int(r[2])))
 with gzip.open(rmsk,'rt') as f:
  for r in csv.reader(f,delimiter='\t'):
   if r[11]=='LINE' and r[5] in lines:lines[r[5]].append((int(r[6]),int(r[7])))
 for r in csv.reader(tssp.open(),delimiter='\t'):
  if r[0] in tss:tss[r[0]].append(int(r[1]))
 vals={k:np.zeros(len(owner),np.int64) for k in ['clean_bp','AT_bp','LINE_bp','TSS_n']}
 with pysam.FastaFile(str(fasta)) as fa:
  for ch in masks:
   ix=np.flatnonzero(a['chromosome']==ch);ss=a['start'][ix];b=np.frombuffer(fa.fetch(ch).upper().encode(),np.uint8);L=len(b);good=np.isin(b,[65,67,71,84])
   for l,h in masks[ch]:good[l:h]=False
   ln=np.zeros(L,bool)
   for l,h in lines[ch]:ln[l:h]=True
   n=L//30000
   for k,x in [('clean_bp',good),('AT_bp',((b==65)|(b==84))&good),('LINE_bp',ln&good)]:
    cnt=x[:n*30000].reshape(n,30000).sum(1);cnt=np.r_[cnt,x[n*30000:].sum()] if L%30000 else cnt;vals[k][ix]=cnt[ss//30000]
   ts=np.array(tss[ch]);ts=ts[good[ts]];vals['TSS_n'][ix]=np.bincount(ts//30000,minlength=(L+29999)//30000)[ss//30000]
 np.testing.assert_array_equal(vals['clean_bp'],a['clean_bp'])
 for k,v in vals.items():seg[k]=np.bincount(owner,weights=v).astype(np.int64)
 for k,num,mul in [('scored_kb','clean_bp',.001),('AT_percent','AT_bp',100),('LINE_percent','LINE_bp',100),('TSS_per_Mb','TSS_n',1e6)]:seg[k]=seg[num]*mul/(seg.clean_bp if k!='scored_kb' else 1)
 seg.to_csv(output/'LAD_context.tsv',sep='\t',index=False);np.savez_compressed(output/'context_bin_features.npz',**vals);print('LAD context complete',flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['signals','decoded','fasta','blacklist','repeatmasker','tss','output']:
  p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 context(a.signals,a.decoded,a.fasta,a.blacklist,a.repeatmasker,a.tss,a.output)
