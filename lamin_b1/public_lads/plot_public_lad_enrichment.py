"""Plot observed and expected activity-gene overlaps with public LAD maps."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt, font_manager
from matplotlib.backends.backend_pdf import PdfPages

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--input',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--fonts',type=Path,required=True)
args=p.parse_args();table=args.input;F=args.output;F.mkdir(parents=True,exist_ok=True)
for fp in args.fonts.glob('arial*.ttf'):font_manager.fontManager.addfont(str(fp))
plt.rcParams.update({'font.family':'Arial','font.size':9,'axes.labelsize':9,'axes.titlesize':9,'xtick.labelsize':8,'ytick.labelsize':9,'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.7})
data=pd.read_csv(table,sep='\t')
labels=['E13 Dorsal NPC','E13 Ventral NPC','30wk Hippocampus','E13 Liver','3T3 MEF','E13 Lung LaminB']
colors=['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']
pdf=F/'public_lad_enrichment.pdf';outputs=[];plotrows=[]
with PdfPages(pdf,metadata={'Title':'Public LAD response-gene permutations: activity-gene enrichment','Author':'','Subject':'Observed and expected activity-gene overlap'}) as book:
    for s,n in [('dPRG',116),('SRG',38),('rPRG',19)]:
        rows=data[data.gene_set==s].set_index('condition').loc[labels]
        assert (rows.genes==n).all()
        obs=rows.observed_percent.to_numpy();mu=rows.null_mean.to_numpy();sd=rows.null_sd.to_numpy()
        fig,ax=plt.subplots(figsize=(180/25.4,99/25.4))
        fig.subplots_adjust(left=.115,right=.982,bottom=.25,top=.885)
        ax.bar(np.arange(6),obs,width=.4 if s=='dPRG' else .6,color=colors,edgecolor='black',linewidth=.9,label='Observed overlap')
        ax.errorbar(np.arange(6),mu,yerr=2*sd,fmt='o',color='black',ms=3.5,elinewidth=1,capsize=3,label='Expected overlap ± 2 SD')
        top=max(max(obs),max(mu+2*sd));ax.set_ylim(0,top+16)
        for i,(_,row) in enumerate(rows.iterrows()):
            ax.text(i,max(obs[i],mu[i]+2*sd[i])+1.5,f'$P$={row.P_unadjusted:.4f}',ha='center',va='bottom',fontsize=8)
            plotrows.append(dict(gene_set=s,map=labels[i],genes=n,observed_percent=obs[i],null_mean=mu[i],null_lower_2SD=mu[i]-2*sd[i],null_upper_2SD=mu[i]+2*sd[i],P_unadjusted=row.P_unadjusted))
        ax.set_ylabel('Gene bodies in LADs (%)')
        ax.set_title(f'{s}s overlap with lamina-associated domains (≥75%; n={n})',pad=8,fontweight='normal')
        ax.set_xticks(np.arange(6),labels,rotation=20,ha='right',rotation_mode='anchor')
        ax.tick_params(width=.7,length=3,pad=3);ax.grid(False)
        ax.legend(loc='upper right',frameon=True,borderpad=.3,handlelength=1.5,labelspacing=.35)
        png=F/f'public_lad_enrichment_{s}.png'
        fig.savefig(png,dpi=220);book.savefig(fig);plt.close(fig);outputs.append(png)
outputs.append(pdf)
values=F/'public_LAD_plot_values.tsv';pd.DataFrame(plotrows).to_csv(values,sep='\t',index=False)
