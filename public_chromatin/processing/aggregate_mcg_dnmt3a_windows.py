"""Aggregate culture/public mCG on the original DNMT3A grid through aligned chain blocks."""
from pathlib import Path
import argparse,gzip,collections,json
from io import StringIO
import numpy as np,pandas as pd
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs',type=Path,required=True)
parser.add_argument('--culture-reports',type=Path,required=True)
parser.add_argument('--public-cache',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--reference-mask',type=Path)
parser.add_argument('--chain',type=Path)
parser.add_argument('--coverage',type=Path)
args=parser.parse_args();inputs=args.inputs;WG=args.culture_reports;public_cache=args.public_cache;output=args.output;output.mkdir(parents=True,exist_ok=True)
chroms=['chr'+str(i) for i in range(1,20)]
grid=pd.read_csv(inputs/'mm9_5kb_windows.tsv',sep='\t');n=len(grid)
off={c:int(np.flatnonzero(grid.chromosome.eq(c))[0]) for c in chroms}
nwin={c:int(grid.chromosome.eq(c).sum()) for c in chroms}
mask=np.load(args.reference_mask or inputs/'reference_mask.npz')
sizes=dict(zip(mask['chromosomes'],mask['lengths']));masks={}
for c in chroms:
    merged=[]
    for a,b in sorted((int(a),int(b)) for typ in ['mask','nonACGT'] for a,b in zip(mask[c+'_'+typ+'_start'],mask[c+'_'+typ+'_end'])):
        if merged and a<=merged[-1][1]:merged[-1][1]=max(b,merged[-1][1])
        else:merged.append([a,b])
    masks[c]=np.asarray(merged,dtype=np.int64)

# Invert aligned chain blocks. Keep only mm10 bases with one mm9 destination.
blocks=collections.defaultdict(list)
chain=args.chain or inputs/'mm9ToMm10.over.chain.gz'
for line in gzip.open(chain,'rt'):
    a=line.split()
    if not a:continue
    if a[0]=='chain':h=a;tp=int(a[5]);qp=int(a[10]);continue
    width=int(a[0])
    if h[2] in chroms and h[7] in chroms:
        assert h[4]=='+'
        rev=h[9]=='-';qs=int(h[8]);lo=qs-qp-width if rev else qp
        blocks[h[7]].append((lo,lo+width,chroms.index(h[2]),tp,int(rev)))
    tp+=width;qp+=width
    if len(a)==3:tp+=int(a[1]);qp+=int(a[2])
mapping={}
for c,rows in blocks.items():
    b=np.asarray(rows,dtype=np.int64);points,idx=np.unique(np.r_[b[:,0],b[:,1]],return_inverse=True)
    count=np.zeros(len(points),dtype=np.int64);identity=count.copy();k=len(b)
    np.add.at(count,idx[:k],1);np.add.at(count,idx[k:],-1)
    np.add.at(identity,idx[:k],np.arange(k)+1);np.add.at(identity,idx[k:],-np.arange(k)-1)
    count=np.cumsum(count);identity=np.cumsum(identity);keep=count[:-1]==1
    segments=np.column_stack([points[:-1][keep],points[1:][keep],identity[:-1][keep]-1])
    mapping[c]=(b,segments)
    assert (segments[:,1]>segments[:,0]).all()
    assert (segments[1:,0]>=segments[:-1,1]).all()

def retain_mm10(c,pos):
    a=masks[c];ix=np.searchsorted(a[:,0],pos,side='right')-1
    return (pos>=0)&(pos//5000<int(sizes[c])//5000)&((ix<0)|(pos>=a[np.maximum(ix,0),1]))

def map_points(c,pos):
    b,seg=mapping[c];idx=np.searchsorted(seg[:,0],pos,side='right')-1;k=np.maximum(idx,0)
    good=(idx>=0)&(pos<seg[k,1]);p=pos[good];bb=b[seg[k[good],2]]
    if not len(p):return np.array([],dtype=np.int64),np.array([],dtype=np.int64)
    pp=bb[:,3]+np.where(bb[:,4]>0,bb[:,1]-1-p,p-bb[:,0])
    # Independent forward round-trip in the original aligned block.
    back=np.where(bb[:,4]>0,bb[:,1]-1-(pp-bb[:,3]),bb[:,0]+(pp-bb[:,3]))
    np.testing.assert_array_equal(back,p)
    codes=bb[:,2];full=pp//5000<np.array([nwin[chroms[int(x)]] for x in np.unique(codes)])[np.searchsorted(np.unique(codes),codes)]
    original_indices=np.flatnonzero(good)[full]
    dest=np.asarray([off[chroms[int(x)]] for x in np.unique(codes)])[np.searchsorted(np.unique(codes),codes)] + pp//5000
    return original_indices,dest[full]

qc={};arrays={}
for label in ['A0','B0','C0','1wk','2wk']:
    cache=output/(label+'_mm9_windows.npz');report=output/(label+'_mapping.json')
    if cache.exists() and report.exists():
        arrays[label]=dict(np.load(cache));qc[label]=json.loads(report.read_text());continue
    met=np.zeros(n,dtype=np.int64);dep=met.copy();sites=met.copy();stats=collections.Counter()
    if label in ['A0','B0','C0']:
        path=next(WG.glob(label+'_*.txt'))
        cutoff=json.loads(((args.coverage or inputs)/f'{label}_coverage.json').read_text())['upper_coverage_exclusive']
        chunks=pd.read_csv(path,sep='\t',header=None,usecols=[0,1,3,4],names=['chr','pos','meth','unmeth'],dtype={'chr':'category','pos':'int64','meth':'int64','unmeth':'int64'},chunksize=1000000)
        def records():
            for chunk in chunks:
                stats['source_records']+=len(chunk)
                chunk['depth']=chunk.meth+chunk.unmeth
                chunk=chunk[(chunk.depth>=4)&(chunk.depth<cutoff)]
                for c,g in chunk.groupby('chr',observed=True):
                    if c in mapping:yield c,g.pos.to_numpy()-1,g.meth.to_numpy(),g.depth.to_numpy()
    else:
        cutoff=json.loads((public_cache/f'public_{label}_summary.json').read_text())['upper_coverage_exclusive']
        def records():
            for c in chroms:
                z=np.load(public_cache/f'public_{label}_CpG_mm10/{c}.npz');stats['cached_records']+=len(z['pos0'])
                take=(z['depth']>=4)&(z['depth']<cutoff)
                yield c,z['pos0'][take].astype(np.int64),z['meth'][take],z['depth'][take]
    for c,pos,mm,nn in records():
        keep=retain_mm10(c,pos);pos,mm,nn=pos[keep],mm[keep],nn[keep]
        stats['eligible_mm10_records']+=len(pos)
        take,ids=map_points(c,pos);mm=mm[take];nn=nn[take]
        stats['retained_mm9_records']+=len(ids);stats['retained_methylated_reads']+=int(mm.sum());stats['retained_total_reads']+=int(nn.sum())
        met+=np.bincount(ids,weights=mm,minlength=n).astype(np.int64)
        dep+=np.bincount(ids,weights=nn,minlength=n).astype(np.int64)
        sites+=np.bincount(ids,minlength=n)
    assert int(sites.sum())==stats['retained_mm9_records']
    assert int(met.sum())==stats['retained_methylated_reads'];assert int(dep.sum())==stats['retained_total_reads']
    assert (met<=dep).all();stats['upper_coverage_exclusive']=cutoff
    stats['excluded_coordinate_or_partial_window']=stats['eligible_mm10_records']-stats['retained_mm9_records']
    arrays[label]={'meth':met,'depth':dep,'sites':sites};np.savez_compressed(cache,**arrays[label])
    qc[label]=dict(stats);report.write_text(json.dumps(qc[label],indent=2)+'\n')
    print(label,'complete',stats['retained_mm9_records'],'CpG records',flush=True)

