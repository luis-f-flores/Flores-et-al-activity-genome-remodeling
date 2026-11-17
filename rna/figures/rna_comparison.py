"""Plot the three published KCl–bicuculline RNA comparisons."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors,font_manager
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--tables',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--fonts',type=Path,required=True)
args=parser.parse_args();T=args.tables;args.output.mkdir(parents=True,exist_ok=True)
d=pd.read_csv(T/'gene_expression_changes.tsv.gz',sep='\t',float_precision='round_trip')
s=pd.read_csv(T/'correlations.tsv',sep='\t',float_precision='round_trip')
for f in args.fonts.glob('arial*.ttf'):font_manager.fontManager.addfont(f)
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'xtick.labelsize':7.5,'ytick.labelsize':7.5})
mm=25.4
fig=plt.figure(figsize=(170/mm,54/mm))
cmap=colors.LinearSegmentedColormap.from_list('gene_density',plt.cm.Blues(np.linspace(.15,1,256)))
for i,tp in enumerate(['1 h','6 h','Recovery']):
 ax=fig.add_axes([(14+i*49)/170,11/54,36/170,36/54])
 z=d[d.timepoint==tp];r=s[(s.timepoint==tp)].iloc[0]
 h=ax.hexbin(z.KCl_TPM_change,z.Bic_TPM_change,gridsize=42,mincnt=1,norm=colors.LogNorm(vmin=1,vmax=4000),cmap=cmap,linewidths=0,extent=(-7.25,7.25,-7.25,7.25),rasterized=True)
 assert int(h.get_array().sum())==len(z)
 ax.axhline(0,color='.65',lw=.5,zorder=0);ax.axvline(0,color='.65',lw=.5,zorder=0)
 ax.set(xlim=(-7.25,7.25),ylim=(-7.25,7.25),xticks=[-5,0,5],yticks=[-5,0,5]);ax.set_aspect('equal')
 ax.set_title(tp,pad=4)
 ax.text(.035,.97,f'ρ = {r.spearman_rho:.3f}\nn = {len(z):,}',transform=ax.transAxes,va='top',fontsize=7.5,bbox={'facecolor':'white','edgecolor':'none','pad':.3})
 if i==0:ax.set_ylabel('6 h Bic / TTX\nlog$_2$ TPM ratio',labelpad=1)
 else:ax.set_yticklabels([])
 ax.tick_params(pad=2)
 fig.text((2+i*49)/170,50/54,chr(ord('G')+i),fontweight='bold',fontsize=11,va='center')
cax=fig.add_axes([156/170,11/54,2/170,36/54]);cb=fig.colorbar(h,cax=cax);cb.ax.set_title('Genes',fontsize=8,pad=4);cb.ax.tick_params(labelsize=7,length=2,pad=2)
fig.savefig(args.output/'Fig_S13_G-I.pdf',dpi=300,facecolor='white')
fig.savefig(args.output/'Fig_S13_G-I.png',dpi=300,facecolor='white')
plt.close(fig)
