"""Bdnf/Mpp7 subtraction tracks on the original 10-kb grid (Figure S8D-E)."""
import os
os.environ.update(MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import numpy as np,pandas as pd,matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from pathlib import Path
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'locus_tracks')
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--fonts',type=Path,default=None)
args=ap.parse_args(); IN=args.inputs; OUT=args.output; OUT.mkdir(parents=True,exist_ok=True)

font_dir=args.fonts
if font_dir is None:
    font_dir=next((p for p in [Path('C:/Windows/Fonts'),Path('/mnt/c/Windows/Fonts'),Path('/usr/share/fonts/truetype/msttcorefonts')] if (p/'arial.ttf').exists()),None)
if font_dir is None: raise SystemExit('Set --fonts to the directory containing Arial .ttf files.')

a=dict(np.load(IN/'supplement.npz'));fine=np.load(IN/'raw_10kb.npz')
ch=fine['chromosome'];st=fine['start'];raw=fine['condition_depth'];groups=[[0,1],[2,3,4],[5,6,7]]
lookup={(str(c),int(s)):i for i,(c,s) in enumerate(zip(ch,st))}
genes=pd.read_csv(IN/'supplement_genes.tsv',sep='	').set_index('gene_name')
views={};arrays={};rows=[]
for gene in ['Bdnf','Mpp7']:
    r=genes.loc[gene];flank=50000 if gene=='Bdnf' else 100000;lo=int(r.start0)-flank;hi=int(r.end0)+flank;x=a[gene+'__B1_starts'];curves={}
    for width in [3,9]:
        half=width//2;cur=np.full((len(x),3),np.nan)
        for j,s in enumerate(x):
            ix=[lookup.get((str(r.chromosome),int(s)+10000*k)) for k in range(-half,half+1)]
            if any(k is None for k in ix):continue
            ip=np.log2(np.mean(2**raw[ix],axis=0));cur[j]=[np.mean(ip[g]) for g in groups]
        curves[width]=cur
    delta=curves[3][:,1:]-curves[3][:,[0]]
    views[gene]=dict(r=r,lo=lo,hi=hi,x=x,delta=delta,tss=int(r.start0) if gene=='Bdnf' else int(r.end0)-1)
    arrays.update({gene+'_starts':x,gene+'_condition_means30':curves[3],gene+'_delta':delta})
    rows.append(dict(gene=gene,chromosome=r.chromosome,lo=lo,hi=hi,gene_start0=int(r.start0),gene_end0=int(r.end0),supported30=int(np.isfinite(curves[3][(x<hi)&(x+10000>lo)]).all(1).sum())))
vals=np.concatenate([d['delta'][(d['x']<d['hi'])&(d['x']+10000>d['lo'])].ravel() for d in views.values()]);ymin=min(-.1,float(np.floor(np.nanmin(vals)*10)/10));ymax=max(.1,float(np.ceil(np.nanmax(vals)*10)/10));print('Subtraction limits',ymin,ymax)

for fn in ['arial.ttf','arialbd.ttf','ariali.ttf']:fm.fontManager.addfont(str(font_dir/fn))
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,'axes.unicode_minus':False})
PW=510.236;PH=748.346;fig=plt.figure(figsize=(PW/72,PH/72));COL=['#0072B2','#D55E00','#4D4D4D'];LS=['-','--','-.']
for i,gene in enumerate(['Bdnf','Mpp7']):
    d=views[gene];left=62.362+i*243.780;wid=184.252
    ax=fig.add_axes([left/PW,213.276/PH,wid/PW,48.189/PH]);ax.axhline(0,color=COL[0],lw=.7)
    ax.axvline(d['tss']/1e6,color='#89949e',lw=.5,ls=':')
    ax.axvspan(d['r'].start0/1e6,d['r'].end0/1e6,color='#f1f1f1',zorder=-3,lw=0)
    for k in range(2):
        ax.plot((d['x']+5000)/1e6,d['delta'][:,k],color=COL[k+1],ls=LS[k+1],lw=.85)
        ax.fill_between((d['x']+5000)/1e6,0,d['delta'][:,k],color=COL[k+1],alpha=.09,lw=0)
    ax.set_xlim(d['lo']/1e6,d['hi']/1e6);ax.set_ylim(ymin,ymax);ax.set_yticks(np.arange(np.ceil(ymin/.4)*.4,ymax+.001,.4));ax.set_ylabel('Δlog2(B1/Input)',fontsize=9,labelpad=3);ax.tick_params(axis='x',bottom=False,labelbottom=False);ax.tick_params(axis='y',length=2,pad=2)
for ax in fig.axes:ax.patch.set_alpha(0)
fig.canvas.draw()

fig.savefig(OUT/'Fig_S08DE.pdf');fig.savefig(OUT/'Fig_S08DE.png',dpi=220);plt.close(fig)
np.savez_compressed(OUT/'subtraction_tracks.npz',**arrays)
pd.DataFrame(rows).to_csv(OUT/'locus_summary.tsv',sep='	',index=False)
