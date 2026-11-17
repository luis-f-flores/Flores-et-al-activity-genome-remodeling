'Draw the promoter/activity-gene timecourse heatmap from the supplied values.'
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt,font_manager
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--input',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'promoter_overlap/activity_gene_timecourse.tsv')
p.add_argument('--output',type=Path,required=True);p.add_argument('--fonts',type=Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(a.input,sep='\t')
for name in ['arial.ttf','arialbd.ttf']:font_manager.fontManager.addfont(str(a.fonts/name))
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.linewidth':.65,'pdf.fonttype':42})
width,height=170,195
fig=plt.figure(figsize=(width/25.4,height/25.4))
def ax(x,y,w,h):return fig.add_axes([x/width,1-(y+h)/height,w/width,h/height])
def txt(x,y,s,**kw):kw.setdefault('va','top');return fig.text(x/width,1-y/height,s,**kw)
classes=['rPRG','dPRG','SRG'];sizes=[19,116,38]
assert list(d.groupby('gene_class',sort=False).size())==sizes and d.gene.nunique()==173
for cl in classes:
 z=d[d.gene_class==cl];assert list(z.gene)==list(z.sort_values(['mCG_0hr_mean','gene'],kind='stable').gene)
mcols=['mCG_0hr_mean']+['mCG_'+t+g+'_mean' for g in ['M','K'] for t in ['1','6','24']]
dcols=['mCG_'+t+'K_minus_'+t+'M' for t in ['1','6','24']]
rcols=['RNA_'+t+'h_log2FC' for t in ['1','6','24']]
for t in ['1','6','24']:
 np.testing.assert_allclose(d['mCG_'+t+'K_minus_'+t+'M'],d['mCG_'+t+'K_mean']-d['mCG_'+t+'M_mean'],atol=1e-12,rtol=0)
assert np.isfinite(d[mcols+dcols+rcols]).all().all()
assert d[dcols].abs().max().max()<25 and d[rcols].abs().max().max()<8
specs=[dict(key='mCG',x=87,w=35,cols=mcols,lo=0,hi=100,cm='YlOrBr'),dict(key='delta',x=128,w=15,cols=dcols,lo=-25,hi=25,cm='RdBu_r'),dict(key='RNA',x=149,w=15,cols=rcols,lo=-8,hi=8,cm='RdBu_r')]
txt(104.5,80,'mCG',ha='center',fontsize=8.5);txt(135.5,80,'ΔmCG',ha='center',fontsize=8.5);txt(156.5,80,'ΔRNA',ha='center',fontsize=8.5)
txt(99.5,86,'Mannitol',ha='center',fontsize=7.7);txt(114.5,86,'KCl',ha='center',fontsize=7.7)
for sp in specs:
 labs=['0','1','6','Recovery','1','6','Recovery'] if sp['key']=='mCG' else ['1','6','Recovery']
 for i,lab in enumerate(labs):txt(sp['x']+(i+.5)*sp['w']/len(labs),101,lab,ha='center',va='bottom',fontsize=7.5,rotation=90 if lab=='Recovery' else 0)
txt(166.5,97,'h',ha='left',fontsize=7.5)
y=103;ims={};display=[]
for cl,n in zip(classes,sizes):
 z=d[d.gene_class==cl];h=n*.34
 txt(84.7,y+h/2,cl+'\nn='+str(n),ha='right',va='center',fontsize=8)
 for sp in specs:
  aa=ax(sp['x'],y,sp['w'],h);mat=z[sp['cols']].to_numpy();cm=plt.get_cmap(sp['cm']).copy();cm.set_bad('#D0D0D0')
  ims[sp['key']]=aa.imshow(mat,aspect='auto',interpolation='nearest',cmap=cm,vmin=sp['lo'],vmax=sp['hi']);aa.set_axis_off()
  for edge in np.arange(.5,mat.shape[1]-.5):aa.axvline(edge,color='white',lw=.65 if sp['key']=='mCG' and edge in [.5,3.5] else .3)
  display.append(dict(gene_class=cl,genes=list(z.gene),measurement=sp['key'],columns=sp['cols'],values=mat.tolist(),limits=[sp['lo'],sp['hi']]))
 y+=h+4
for key,x,w,unit,ticks in [('mCG',85,23,'%', [0,50,100]),('delta',114,23,'ΔmCG (%)',[-25,0,25]),('RNA',143,23,'log2 fold change',[-8,0,8])]:
 cb=fig.colorbar(ims[key],cax=ax(x,179,w,2),orientation='horizontal',ticks=ticks);cb.outline.set_linewidth(.5);cb.ax.tick_params(labelsize=7.5,pad=1.8,length=1.8);txt(x+w/2,187,unit,ha='center',fontsize=7.2)

fig.savefig(a.output/'promoter_activity_timecourse.pdf',metadata={'CreationDate':None,'ModDate':None})
fig.savefig(a.output/'promoter_activity_timecourse.png',dpi=200)
plt.close(fig)
print('Reproduced promoter/activity-gene timecourse component.')
