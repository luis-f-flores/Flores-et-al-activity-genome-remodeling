"""Draw Figure 1D from condition-mean B1 measurements."""
import os
os.environ.update(MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import argparse,json,csv
import numpy as np
from matplotlib import pyplot as plt,font_manager

p=argparse.ArgumentParser()
p.add_argument('--out',type=Path,required=True)
p.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'b1_spread')
p.add_argument('--font-dir',type=Path,default=Path('C:/Windows/Fonts'))
args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
source=args.inputs/'D_bin_means.npz'
m=np.load(source)['means'];assert m.ndim==2 and m.shape[1]==3 and np.isfinite(m).all()
sd=m.std(axis=0,ddof=1)
np.testing.assert_allclose(sd,np.sqrt(((m-m.mean(0))**2).sum(0)/(len(m)-1)),atol=1e-14,rtol=0)
pairs=[(0,1),(0,2),(1,2)];names=['TTX','Bic','Bic+Nimo']
rows=[{'reference':names[i],'target':names[j],'n_bins':len(m),'SD_reference':float(sd[i]),'SD_target':float(sd[j]),'SD_ratio':float(sd[j]/sd[i])} for i,j in pairs]
for name in ['arial.ttf','arialbd.ttf']:
    if (args.font_dir/name).is_file():font_manager.fontManager.addfont(str(args.font_dir/name))
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.labelsize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'savefig.facecolor':'white'})
width,height=510.235992431641,115
fig=plt.figure(figsize=(width/72,height/72));labels=['TTX','Bicuculline','Bic + nimodipine']
for k,((i,j),row) in enumerate(zip(pairs,rows)):
    ax=fig.add_axes([(73+158.74*k)/width,34/height,74/width,74/height])
    ax.hexbin(m[:,i],m[:,j],gridsize=65,mincnt=1,bins='log',cmap='Greys',linewidths=0,rasterized=True,extent=(-2,3.2,-2,3.2))
    ax.plot([-2,3.2],[-2,3.2],color='#0072B2',ls='--',lw=.8)
    ax.text(.035,.955,'SD ratio = %.2f'%row['SD_ratio'],transform=ax.transAxes,ha='left',va='top',fontsize=8.5)
    ax.set(xlim=(-2,3.2),ylim=(-2,3.2),xticks=[-2,0,2],yticks=[-2,0,2],xlabel=labels[i]+'\nlog2(B1/Input)',ylabel=labels[j]+'\nlog2(B1/Input)')
    ax.tick_params(length=2.5,pad=2);ax.xaxis.labelpad=2;ax.yaxis.labelpad=2;ax.set_aspect('equal',adjustable='box')
for ext in ['svg','pdf','png']:fig.savefig(args.out/('Fig_1D.'+ext),dpi=450)
plt.close(fig)
with (args.out/'SD_ratios.tsv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
(args.out/'RESULTS.json').write_text(json.dumps({'ratios':rows,'definition':'SD(y-axis condition)/SD(x-axis condition)','descriptive':True},indent=2)+'\n',encoding='utf-8')
print('SD ratios:',[round(row['SD_ratio'],2) for row in rows])
