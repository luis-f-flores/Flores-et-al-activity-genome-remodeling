"""Ccbe1 and Grin2a B1/RNA locus profiles (Figure 3C-D)."""
import os
os.environ.update(MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
import numpy as np,pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator,FormatStrFormatter
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

for fn in ['arial.ttf','arialbd.ttf','ariali.ttf']:fm.fontManager.addfont(str(font_dir/fn))
v=pd.read_csv(IN/'main_genes.tsv',sep='	').set_index('gene_name')
gm=pd.read_csv(IN/'gene_metadata.tsv',sep='	');ex=pd.read_csv(IN/'exons.bed',sep='	',header=None)
st=pd.read_csv(IN/'main_summary.tsv',sep='	').set_index('gene')
a=np.load(IN/'main.npz');win=np.load(IN/'window_bins.npy')*10
edges=np.empty(len(win)+1);edges[1:-1]=np.sqrt(win[:-1]*win[1:]);edges[0]=win[0]**2/edges[1];edges[-1]=win[-1]**2/edges[-2]
view={}
for gene in ['Ccbe1','Grin2a']:
 r=v.loc[gene];q=st.loc[gene]
 loci=set(gm[(gm.gene_name==gene)&(gm.chrom==r.chromosome)&(gm.start0<r.end0)&(gm.end0>r.start0)].locus_id)
 view[gene]=dict(r=r,lo=int(q.view_start0),hi=int(q.view_end0),x=a[gene+'_starts'],delta=a[gene+'_delta'],z=a[gene+'_domainograms'],rx=a[gene+'_RNA_starts'],rna=a[gene+'_RNA_CPM'],rna_mask=a[gene+'_RNA_display_mask'],rmax=q.RNA_ymax,rlo=int(q.RNA_start0),rhi=int(q.RNA_end0),hmm_starts=a[gene+'_HMM_starts'],hmm_states=a[gene+'_HMM_states'],ex=ex[ex[3].isin(loci)])
plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.55,'axes.unicode_minus':False})
MM=72/25.4;PW=510.236/MM;PH=313/MM;fig=plt.figure(figsize=(510.236/72,313/72));COL=['#0072B2','#D55E00','#4D4D4D'];LS=['-','--','-.'];SC=['#E1E4E9','#71BBDD','#164A63']
cm=LinearSegmentedColormap.from_list('B1_change',['#762A83','white','#1B7837'],N=257);cm.set_bad('#D9D9D9')
def text(x,y,s,fs=8,**kw):return fig.text(x/PW,1-y/PH,s,fontsize=fs,va='baseline',zorder=1000,clip_on=False,**kw)
def ax(x,y,w,h):return fig.add_axes([x/PW,1-(y+h)/PH,w/PW,h/PH])
def line(x,y,w,c,ls='-'):
    fig.add_artist(Line2D([x/PW,(x+w)/PW],[1-y/PH]*2,color=c,ls=ls,lw=.8))
def square(x,y,c):fig.add_artist(Rectangle((x/PW,1-(y+2)/PH),2/PW,2/PH,facecolor=c,edgecolor='none'))
for x,c,ls,label in [(16,COL[0],LS[0],'TTX'),(48,COL[1],LS[1],'Bicuculline'),(92,COL[2],LS[2],'Bicuculline + nimodipine')]:line(x,10,4,c,ls);text(x+6,11,label,8.5)
text(15,16,'Domainogram',8);square(39,13.5,'#762A83');text(43,16,'Lower B1',8);square(61,13.5,'#1B7837');text(65,16,'Higher B1',8)
text(95,16,'HMM',8);square(105,13.5,SC[0]);text(109,16,'Non-B1',8);square(128,13.5,SC[1]);text(132,16,'Intermediate',8);square(160,13.5,SC[2]);text(164,16,'High',8)
for gi,gene in enumerate(['Ccbe1','Grin2a']):
    d=view[gene];r=d['r'];left=15+gi*94;wid=65;axes=[]
    text(3+gi*94,3.5,'CD'[gi],11,fontweight='bold');text(8+gi*94,3.5,gene,9,fontstyle='italic')
    for ci,y in [(0,23),(1,36)]:
        text(left,y-1.5,'Bicuculline - TTX' if ci==0 else 'Bicuculline + nimodipine - TTX',8)
        aa=ax(left,y,wid,8);axes.append(aa);aa.pcolormesh(np.r_[d['x'],d['x'][-1]+10000]/1e6,edges,d['z'][ci],cmap=cm,vmin=-1,vmax=1,rasterized=True);aa.set_yscale('log');aa.minorticks_off();aa.set_yticks([30,990],['30','990']);aa.set_ylabel('Window\n(kb)',fontsize=8,labelpad=0);aa.tick_params(axis='y',length=0,pad=2);aa.spines['bottom'].set_visible(False)
    aa=ax(left,47,wid,19);axes.append(aa);aa.axhline(0,color=COL[0],lw=.75,zorder=0)
    for k in range(2):aa.plot((d['x']+5000)/1e6,d['delta'][:,k],color=COL[k+1],ls=LS[k+1],lw=.9);aa.fill_between((d['x']+5000)/1e6,0,d['delta'][:,k],color=COL[k+1],alpha=.09,lw=0)
    aa.axvspan(r.start0/1e6,r.end0/1e6,color='#f1f1f1',lw=0,zorder=-2);aa.set_ylim(-.5,.9);aa.set_yticks([-.5,0,.3,.6,.9]);aa.set_ylabel('Δlog2(B1/Input)',fontsize=8.5,labelpad=2);aa.spines['bottom'].set_visible(False);aa.tick_params(axis='y',length=2,pad=2)
    assert np.nanmin(d['delta'])>=-.5 and np.nanmax(d['delta'])<=.9
    aa.text(.99,.97,'Change from TTX',ha='right',va='top',fontsize=7.5,transform=aa.transAxes)
    aa=ax(left,68,wid,3.6);axes.append(aa);hmm_ax=aa
    for k in range(3):
        for j,s0 in enumerate(d['hmm_starts']):
            aa.add_patch(Rectangle((s0/1e6,2-k-.4),.03,.8,facecolor=SC[d['hmm_states'][j,k]],edgecolor='none'))
        fig.add_artist(Rectangle(((left-1.9)/PW,1-(68+k*1.2+.95)/PH),1.1/PW,.95/PH,facecolor=COL[k],edgecolor='none'))
    aa.set_ylim(-.5,2.5);aa.set_yticks([]);aa.spines['left'].set_visible(False);aa.spines['bottom'].set_visible(False);text(left-2.5,70.5,'HMM',8,ha='right')
    text(left,76.5,'RNA (CPM)',8.5)
    text(left,80,f'log2FC {r.fc_Bic_vs_TTX:+.2f}; q={r.q_Bic_vs_TTX:.3g}',7.8)
    rna_axes=[]
    for k in range(3):
        aa=ax(left,82+k*5.5,wid,5);rna_axes.append(aa)
        yy=np.where(d['rna_mask'],d['rna'][:,k],np.nan)
        np.testing.assert_array_equal(yy[d['rna_mask']],d['rna'][d['rna_mask'],k])
        aa.fill_between(d['rx']/1e6,0,yy,step='post',color=COL[k],lw=0)
        aa.set_ylim(0,d['rmax']);aa.set_yticks([d['rmax']],[f"{d['rmax']:.1f}"])
        aa.tick_params(axis='y',length=0,pad=2);aa.spines['left'].set_visible(False);aa.spines['bottom'].set_color('#8c8c8c');aa.spines['bottom'].set_linewidth(.3)
    aa=ax(left,99,wid,3);rna_axes.append(aa);aa.plot([r.start0/1e6,r.end0/1e6],[0,0],color='#555555',lw=.6)
    for er in d['ex'].itertuples(index=False,name=None):aa.add_patch(Rectangle((int(er[1])/1e6,-.2),(int(er[2])-int(er[1]))/1e6,.4,color='#555555',lw=0))
    dx=(1 if r.strand=='+' else -1)*(d['rhi']-d['rlo'])/1e6*.06;aa.plot([r.tss0/1e6]*2,[.2,.9],color='#555555',lw=.7);aa.annotate('',xy=(r.tss0/1e6+dx,.9),xytext=(r.tss0/1e6,.9),arrowprops=dict(arrowstyle='-|>',color='#555555',lw=.7),annotation_clip=False)
    aa.set_ylim(-.35,1.2);aa.set_yticks([]);aa.spines['left'].set_visible(False);aa.set_xlabel(r.chromosome+' (Mb)',fontsize=8,labelpad=1);ticks=MaxNLocator(4).tick_values(d['rlo']/1e6,d['rhi']/1e6);aa.set_xticks(ticks[(ticks>=d['rlo']/1e6)&(ticks<=d['rhi']/1e6)]);aa.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    for a in axes:
        a.set_xlim(d['lo']/1e6,d['hi']/1e6);a.tick_params(labelsize=8)
        a.tick_params(axis='x',bottom=False,labelbottom=False)
    for a in rna_axes:
        a.set_xlim(d['rlo']/1e6,d['rhi']/1e6);a.tick_params(labelsize=8)
        if a is not aa:a.tick_params(axis='x',bottom=False,labelbottom=False)
for axx in fig.axes:
    axx.patch.set_alpha(0)
    for tt in axx.findobj(plt.Text):tt.set_clip_on(False)
fig.canvas.draw();ren=fig.canvas.get_renderer()
for t in fig.findobj(plt.Text):
    if t.get_visible() and t.get_text():
        b=t.get_window_extent(ren);assert b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+1,(t.get_text(),b.bounds)
fig.savefig(OUT/'Fig_3CD.pdf');fig.savefig(OUT/'Fig_3CD.png',dpi=220);plt.close(fig)
