'Draw promoter-DMG count and locus panels. TPM panels use methylation/figures/promoter_methylation.py.'

from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)

from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle,Rectangle
from scipy.stats import t as student_t
W=Path(gb_path('wgbs_package'))
R=W/'runs/2026-10-14_RNA_DMR_main_supp';D=R/'69_main3';T=W/'runs/2026-09-12_fig4_gene_mapping_ols/tables'
O=Path(gb_path('figures/WGBS/runs/2026-10-14_RNA_DMR_main_supp/figures/2026-10-18_main3_compact'));O.mkdir(parents=True,exist_ok=True)
for fn in ['arial.ttf','arialbd.ttf','ariali.ttf']:font_manager.fontManager.addfont(str(Path(gb_path('@fonts'))/fn))
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':2,'ytick.major.size':2})
PW,PH=180,246;fig=plt.figure(figsize=(PW/25.4,PH/25.4))
def ax(x,y,w,h):return fig.add_axes([x/PW,1-(y+h)/PH,w/PW,h/PH])
def text(x,y,s,size=9,**kw):return fig.text(x/PW,1-y/PH,s,fontsize=size,**kw)
def title(letter,s,x,y):text(x,y,letter,11,fontweight='bold');text(x+6,y,s,9)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:pd.read_csv(p,sep='\t')
full={t:read(T/f'20_DMG_gene_matched_{t}hr.tsv') for t in [1,6,24]}
COL={'0 h':'#626970','Mannitol':'#287BB5','KCl':'#CA474E'};hypo='#278FC4';hyper='#E64B66'
titles=[];replay=[]
title('A','Promoter DMGs',3,6)
labels={1:['Rcan2','Olfr545','Mir692-2','Mir692-3','Dusp12','Scgb2b24','Prkar1a','Mir7043','Olfr1489','Ocln'],6:['9630028H03Rik','Pard3bos2','Tdpoz8','Efemp1','Mir7675','Kti12','Mir692-3','Mir692-2','Mblac1','Speer4e'],24:['Gucy1b1','Pnpt1','F730043M19Rik','Itgb1bp2','Pla2g4e','Fktn','Chadl','Gm5800','Zfp52','Napepld']}
for j,t in enumerate([1,6,24]):
 a=ax(14+59*j,14,44,29);d=full[t];y=-np.log10(d.qvalue);ymax={1:12,6:22,24:12}[t]
 a.scatter(d.meth_diff,y,c=np.where(d.meth_diff<0,hypo,hyper),s=7,edgecolors='none',zorder=2)
 a.axhline(2,color='#888',ls=(0,(2,3)),lw=.5)
 for v in [-25,25]:a.axvline(v,color='#888',ls='--',lw=.5)
 a.axvline(0,color='#aaa',ls=(0,(2,3)),lw=.45)
 a.set(xlim=(-85,85),ylim=(1.5,ymax),xticks=[-80,-40,0,40,80],yticks={1:[3,6,9,12],6:[5,10,15,20],24:[3,6,9,12]}[t])
 a.set_title(('Recovery' if t==24 else f'{t} h'),pad=3);a.tick_params(pad=1)
 if j==0:a.set_ylabel(r'−log$_{10}$(FDR)',labelpad=1)
 for sign in [-1,1]:
  named=d[d.gene_name.isin(labels[t]) & (np.sign(d.meth_diff)==sign)].copy();named['yy']=-np.log10(named.qvalue);named=named.sort_values(['yy','gene_name'])
  targets=np.maximum(named.yy.to_numpy(),np.linspace(2.8,ymax*.7,len(named)))
  min_gap=(ymax-1.5)*.115
  for k in range(1,len(targets)):targets[k]=max(targets[k],targets[k-1]+min_gap)
  if len(targets) and targets[-1]>ymax*.94:targets-=targets[-1]-ymax*.94
  for (_,r),yy in zip(named.iterrows(),targets):
   xtext=-81 if sign<0 else 81
   a.annotate(r.gene_name,(r.meth_diff,r.yy),xytext=(xtext,yy),fontsize=7.2,fontstyle='italic',ha='left' if sign<0 else 'right',va='center',arrowprops=dict(arrowstyle='-',lw=.35,color='#666'),bbox=dict(facecolor='white',alpha=.8,edgecolor='none',pad=.15),zorder=3)
 a.set_xlabel('ΔmCG (percentage points)',labelpad=1,fontsize=8)
 assert len(d)=={1:81,6:87,24:90}[t]

title('B','Promoter DMGs and DEGs',3,55)
for x,c,s in [(119,'#F29B38','DMGs'),(146,'#9472BB','DEGs')]:
 fig.add_artist(Circle((x/PW,1-54/PH),1.15/PW,transform=fig.transFigure,facecolor='none',edgecolor=c,lw=.9));text(x+3,55,s,8)
counts=read(T/'20_fig4B_counts_verified.tsv')
for j,t in enumerate([1,6,24]):
 a=ax(6+59*j,60,50,20);a.set(xlim=(0,50),ylim=(0,25));a.set_aspect('equal');a.axis('off')
 a.add_patch(Circle((12,11.7),7.3,fill=False,ec='#F29B38',lw=.85));a.add_patch(Circle((29,11.7),11.1,fill=False,ec='#9472BB',lw=.85))
 row=counts[counts.timepoint==f'{t}hr'].iloc[0];shared=int(row.corrected_overlap)
 a.text(10.9,11.6,str(int(row.all_DMGs)-shared),ha='center',va='center',fontsize=9)
 a.text(30,11.6,str(int(row.DEGs)-shared),ha='center',va='center',fontsize=9)
 a.annotate(str(shared),xy=(18.5,13),xytext=(18.5,23.0),ha='center',fontsize=8,arrowprops=dict(arrowstyle='-',lw=.5))
 a.text(33,24,('Recovery' if t==24 else f'{t} h'),ha='center',fontsize=9)

title('E','6 h hypo-DMGs · SynGO',94,137)


tracks=read(D/'track_plot_values.tsv');rna=read(D/'locus_RNA_values.tsv');ex=read(W/'rna_seq/gene_layers_mm10/exons.bed') if False else pd.read_csv(W/'rna_seq/gene_layers_mm10/exons.bed',sep='\t',header=None,names=['chr','start','end','gene','score','strand'])
genes=pd.read_csv(W/'rna_seq/gene_layers_mm10/genes.bed',sep='\t',header=None,names=['chr','start','end','gene','score','strand'])
for letter,gene,chrom,lo,hi,left in [('F','Ccdc25','chr14',65835120,65839660,13),('G','Serpinb8','chr1',107588762,107592291,104)]:
 text(left-10,186,letter,11,fontweight='bold');text(left-4,186,gene,9,fontstyle='italic');sub=tracks[tracks.gene==gene];dd=full[6][full[6].gene_name==gene].iloc[0];color='#7F4A91' if dd.meth_diff>0 else '#249A8C';width=49
 for gi,(group,col) in enumerate(COL.items()):
  a=ax(left,195+gi*10,width,5.6);d=sub[sub.group==group]
  a.axvspan(dd.prom_start,dd.prom_end,color=color,alpha=.075,lw=0)
  a.vlines(d.position1,0,d.mCG,color=col,lw=.38);a.vlines(d.position1,-2.5,0,color=col,lw=.3)
  a.set(xlim=(lo,hi),ylim=(-4,100),yticks=[0,100]);a.tick_params(axis='x',bottom=False,labelbottom=False);a.tick_params(axis='y',labelsize=8,pad=1)
  a.text(.01,1.12,group if gi==0 else group+' · 6 h',transform=a.transAxes,color=col,fontsize=8.5,va='bottom')
  if gi==1:a.set_ylabel('mCG (%)',labelpad=1,fontsize=8)
 a=ax(left,224,width,3);a.set(xlim=(lo,hi),ylim=(0,1));a.axis('off');a.add_patch(Rectangle((dd.prom_start,.25),dd.prom_end-dd.prom_start+1,.5,color=color,lw=0));a.text(.01,1.2,'Hyper-DMR' if dd.meth_diff>0 else 'Hypo-DMR',transform=a.transAxes,fontsize=8,color=color)
 a=ax(left,231,width,4);a.set(xlim=(lo,hi),ylim=(0,1),yticks=[]);a.spines['left'].set_visible(False)
 gg=genes[genes.gene==gene].iloc[0];assert gg.strand=='+';tss=gg.start+1
 a.plot([max(tss,lo),hi],[.45,.45],color='#333',lw=.7)
 for e in ex[(ex.gene==gene)&(ex.start<hi)&(ex.end>=lo)].itertuples():a.add_patch(Rectangle((max(e.start+1,lo),.20),min(e.end,hi)-max(e.start+1,lo)+1,.45,color='#333',lw=0))
 a.annotate('',xy=(tss+(hi-lo)*.10,.97),xytext=(tss,.97),arrowprops=dict(arrowstyle='->',lw=.7))
 ticks=np.linspace(lo,hi,3);a.set_xticks(ticks,[f'{v/1e6:.4f}' for v in ticks]);a.tick_params(axis='x',length=2,pad=1,labelsize=7.8);a.set_xlabel(chrom+' (Mb; mm10)',fontsize=8.5,labelpad=1)
 a=ax(left+61,200,12,23);d=rna[rna.gene==gene].sort_values('time_h');a.axhline(0,color='#777',lw=.6);a.bar(range(3),d.log2FC,color=COL['Mannitol'],width=.68)
 ymax=1 if gene=='Ccdc25' else 6;ymin=-1 if gene=='Ccdc25' else -3
 a.set(ylim=(ymin,ymax),yticks=[-1,0,1] if gene=='Ccdc25' else [-2,0,2,4,6]);a.set_xticks(range(3),['1','6','Recovery']);a.tick_params(pad=1,labelsize=8);a.set_xlabel('Time (h)',fontsize=8,labelpad=1);a.set_ylabel(r'RNA log$_2$ FC',fontsize=8,labelpad=1)
 for k,r in enumerate(d.itertuples()):a.text(k,r.log2FC+(.05 if r.log2FC>0 else -.025)*(ymax-ymin),f'{r.padj:.3g}',ha='center',va='bottom' if r.log2FC>0 else 'top',fontsize=7.3)
fig.savefig(D/'main3_native_base.pdf');fig.savefig(D/'main3_native_base.png',dpi=240);plt.close(fig)
