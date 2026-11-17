'Draw DMR PCA, locus, heatmap and annotation panels.'
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'wgbs_figures')
ap.add_argument('--output',type=Path,required=True);ap.add_argument('--fonts',type=Path,required=True)
args=ap.parse_args();M=R=args.inputs;OUT=args.output;OUT.mkdir(parents=True,exist_ok=True)
for n in ['arial.ttf','arialbd.ttf','ariali.ttf']:font_manager.fontManager.addfont(str(args.fonts/n))
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.titlesize':9,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
C={'0':'#626970','M':'#287BB5','K':'#CA474E'}
PW,PH=180,238
fig=plt.figure(figsize=(PW/25.4,PH/25.4))
def ax(x,y,w,h):return fig.add_axes([x/PW,1-(y+h)/PH,w/PW,h/PH])
def text(x,y,t,size=9,**kw):fig.text(x/PW,1-y/PH,t,fontsize=size,**kw)
def title(letter,title_,x,y):text(x,y,letter,11,fontweight='bold');text(x+6,y,title_,9)
scores=pd.read_csv(M/'PCA_scores_PC1_PC2_PC3.tsv',sep='\t');var=pd.read_csv(M/'PCA_variance.tsv',sep='\t').variance_percent.to_numpy()
a=ax(102,62,68,35);markers={0:'o',1:'s',6:'^',24:'D'}
for row in scores.itertuples():
 tr='0' if row.sample.endswith('0') else row.sample[-1];t=0 if tr=='0' else int(row.sample[1:-1])
 a.scatter(row.PC1/1000,row.PC2/1000,s=23,c=C[tr],marker=markers[t],edgecolor='white',linewidth=.35)
a.set(xlabel=f'PC1 ({var[0]:.1f}%)',ylabel=f'PC2 ({var[1]:.1f}%)');a.margins(.13)
a.xaxis.labelpad=1;a.yaxis.labelpad=1;a.tick_params(pad=1,labelsize=8)
text(150,60,'Scores ×10³',8)
handles=[Line2D([],[],marker='o',linestyle='',color=C[t],label=l) for t,l in [('0','0 h'),('M','Mannitol'),('K','KCl')]]
handles += [Line2D([],[],marker=markers[t],linestyle='',color='#444444',label=('Recovery' if t==24 else f'{t} h')) for t in [1,6,24]]
fig.legend(handles=handles,loc='center',bbox_to_anchor=(137/PW,1-109/PH),frameon=False,ncol=3,fontsize=9,handletextpad=.2,columnspacing=.8,labelspacing=.2)
# Three DMR loci with observed-replicate means and shared 0–100% axes.
d=pd.read_csv(R/'three_example_plot_values.tsv',sep='\t');dm=pd.read_csv(R/'three_example_DMRs.tsv',sep='\t')
ex=pd.read_csv(M/'exons.bed',sep='\t',header=None,names=['chr','start','end','gene','score','strand'])
for gene,chrom,lo,hi,t,tss,left in [('Adssl1','chr12',112628751,112630000,1,None,13),('Bdnf','chr2',109692151,109693150,6,109692436,73),('Styx','chr14',45349100,45353609,24,45351186,133)]:
 width=44;sub=d[d.gene==gene];dd=dm[dm.gene==gene];color='#249A8C' if gene=='Bdnf' else '#7F4A91'
 text(left,118,gene,9,fontstyle='italic')
 if gene=='Bdnf':text(left+8.5,118,'exon IV region',9)
 text(left+width,118,('Recovery' if t==24 else f'{t} h'),9,ha='right')
 for gi,(label,key,col) in enumerate([('0 h','0 h',C['0']),(f'Mannitol','Mannitol',C['M']),(f'KCl','KCl',C['K'])]):
  a=ax(left,123+gi*11,width,7)
  for row in dd.itertuples():a.axvspan(row.start,row.end,color=color,alpha=.075,lw=0,zorder=0)
  group=sub[sub.group==key];assert group.mCG.between(0,100).all()
  a.vlines(group.position1,0,group.mCG,color=col,lw=.45)
  a.vlines(group.position1,-2.5,0,color=col,lw=.35)
  a.set(xlim=(lo,hi),ylim=(-4,100),yticks=[0,100]);a.tick_params(axis='x',bottom=False,labelbottom=False);a.tick_params(axis='y',labelsize=8,pad=1,length=2)
  a.text(.015,1.05,label,transform=a.transAxes,color=col,fontsize=9,va='bottom')
  if gi==1:a.set_ylabel('mCG (%)',fontsize=9,labelpad=1)
 a=ax(left,155,width,4);a.set(xlim=(lo,hi),ylim=(0,1));a.axis('off')
 for row in dd.itertuples():a.add_patch(Rectangle((row.start,.2),row.end-row.start+1,.5,color=color,lw=0))
 a.text(.01,-.7,'Hypo-DMR' if gene=='Bdnf' else 'Hyper-DMR',transform=a.transAxes,fontsize=9,color=color)
 a=ax(left,164,width,5);a.set(xlim=(lo,hi),ylim=(0,1),yticks=[]);a.spines['left'].set_visible(False)
 a.plot([max(tss,lo) if tss is not None else lo,hi],[.45,.45],color='#333',lw=.7)
 for e in ex[(ex.gene==gene)&(ex.start<hi)&(ex.end>=lo)].itertuples():a.add_patch(Rectangle((max(e.start+1,lo),.2),min(e.end,hi)-max(e.start+1,lo)+1,.45,color='#333',lw=0))
 if tss is not None:
  a.annotate('',xy=(tss+(hi-lo)*.09,.96),xytext=(tss,.96),arrowprops=dict(arrowstyle='->',lw=.7))
 else:
  # Adssl1 intronic fragment, positive strand; its TSS is outside this window.
  for frac in [.28,.62]:
   a.annotate('',xy=(lo+(hi-lo)*(frac+.09),.45),xytext=(lo+(hi-lo)*frac,.45),arrowprops=dict(arrowstyle='->',lw=.7))
 ticks=np.array([lo,hi]);a.set_xticks(ticks);a.set_xticklabels([f'{v/1e6:.4f}' for v in ticks],fontsize=8)
 a.get_xticklabels()[0].set_ha('left');a.get_xticklabels()[-1].set_ha('right')
 a.tick_params(axis='x',length=2,pad=1);a.set_xlabel(chrom+' (Mb; mm10)',fontsize=9,labelpad=1)
z=np.load(M/'DMR_sample_matrix.npz');ids=list(z['samples']);m=z['mCG_percent'];keep=np.isfinite(m).all(axis=1)&(np.nanstd(m,axis=1)>0);x=m[keep];x=x-x.mean(axis=1,keepdims=True)
ordering=pd.read_csv(M/'heatmap_cluster_order.tsv.gz',sep='\t');assert np.array_equal(ordering.region_id.to_numpy(),z['region_ids'][keep])
order=np.argsort(ordering.display_order.to_numpy());assert np.array_equal(np.sort(order),np.arange(len(x)))
sizes=ordering.groupby('display_group',sort=True).size().tolist()
for left,values,counts,limit,ylabel,legend in [(15,x[order],sizes,20,'229,544 DMR regions','ΔmCG (%)')]:
 a=ax(left,194,68,30);im=a.imshow(values,aspect='auto',interpolation='nearest',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
 a.set(xticks=[],yticks=[]);a.set_ylabel(ylabel,fontsize=9,labelpad=2)
 for cut in [2.5,5.5,8.5,11.5,14.5,17.5]:a.axvline(cut,color='white',lw=.5)
 for cut in np.cumsum(counts)[:-1]:a.axhline(cut-.5,color='white',lw=.6)
 for j,label in enumerate(['0','1','6','Recovery','1','6','Recovery']):a.text((j*3+1.5)/21,1.015,label,transform=a.transAxes,ha='center',fontsize=7.5,rotation=90 if label=='Recovery' else 0)
 for lo,hi,label,col in [(0,3,'0 h',C['0']),(3,12,'Mannitol',C['M']),(12,21,'KCl',C['K'])]:
  a.plot([lo/21,hi/21],[1.42,1.42],transform=a.transAxes,color=col,lw=2,clip_on=False)
  a.text((lo+hi)/42,1.48,label,transform=a.transAxes,ha='center',fontsize=9,color=col)
 cb=fig.colorbar(im,cax=ax(left+12,228,40,1.7),orientation='horizontal',ticks=[-limit,0,limit],extend='both');cb.ax.tick_params(labelsize=8,length=2,pad=1);cb.set_label(legend,fontsize=9,labelpad=1)
fig.savefig(OUT/'Fig_4_DEF.pdf');plt.close(fig)

PW,PH=180,118
fig=plt.figure(figsize=(PW/25.4,PH/25.4))
title('A','Global mCG',4,7)
title('B','DMR annotations',75,7)
pies=json.loads((M/'original_DMR_annotation_pies.json').read_text())
for ri,direction in enumerate(['hypo','hyper']):
 text(125,18+46*ri,'Hypomethylated DMRs' if ri==0 else 'Hypermethylated DMRs',9,ha='center')
 for ci,t in enumerate([1,6,24]):
  row=next(r for r in pies if r['direction']==direction and r['timepoint']==t)
  keys=['promoter_dmrs','gene_body_dmrs','intergenic_dmrs','CGI']
  vals=[row['counts'][k] for k in keys]
  assert sum(vals)==row['displayed_total']==row['sum_of_four_categories']
  assert sum(vals[:3])==row['unique_genomic_DMRs']
  for k,v in zip(keys,vals):assert round(100*v/sum(vals),2)==row['percentages'][k]
  left=75+33*ci;top=25+46*ri
  a=ax(left,top,31,28)
  wedges,_=a.pie(vals,colors=['white','#a578b6','#6f80c1','black'],startangle=90,counterclock=False,wedgeprops={'edgecolor':'black','linewidth':.5},radius=1)
  for j in [1,2]:
   angle=np.deg2rad((wedges[j].theta1+wedges[j].theta2)/2)
   a.text(.58*np.cos(angle),.58*np.sin(angle),f'{100*vals[j]/sum(vals):.2f}%',color='white',ha='center',va='center',fontsize=8)
  for j,xx in [(0,.53),(3,-.53)]:
   angle=np.deg2rad((wedges[j].theta1+wedges[j].theta2)/2)
   a.annotate(f'{100*vals[j]/sum(vals):.2f}%',xy=(.99*np.cos(angle),.99*np.sin(angle)),xytext=(xx,1.22),ha='center',fontsize=8,arrowprops={'arrowstyle':'-','lw':.5,'color':'#444'})
  a.set(xlim=(-1.13,1.13),ylim=(-1.1,1.45));a.axis('off')
  text(left+15.5,top-1,('Recovery' if t==24 else f'{t} h'),8,ha='center')
  text(left+15.5,top+29,f"{sum(vals):,} annotations",8,ha='center')
from matplotlib.patches import Patch
ph=[Patch(facecolor=col,edgecolor='black',linewidth=.5,label=lab) for lab,col in [('Promoter','white'),('Gene body','#a578b6'),('Intergenic','#6f80c1'),('CpG island','black')]]
fig.legend(handles=ph,loc='center',bbox_to_anchor=(125/PW,1-111/PH),ncol=2,frameon=False,fontsize=8,handlelength=1,columnspacing=1.2,handletextpad=.4)
fig.savefig(OUT/'Fig_S12B.pdf');plt.close(fig)
