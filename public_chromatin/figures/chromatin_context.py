"""Plot Figure S10 from the recomputed chromatin-context tables."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--tables',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--fonts',type=Path,required=True)
args=parser.parse_args();T=args.tables;FO=args.output;FO.mkdir(parents=True,exist_ok=True)
inputs={'A':T/'ranked_LAD_groups.tsv','rho':T/'B1_correlations.tsv','B':T/'paired_body_profiles.tsv','C':T/'mean_DNMT3A_rank_groups.tsv','DE_trends':T/'quantile_trends.tsv','DE_correlations':T/'association_summary.tsv'}
g=pd.read_csv(inputs['A'],sep='\t');r=pd.read_csv(inputs['rho'],sep='\t');b=pd.read_csv(inputs['B'],sep='\t');m=pd.read_csv(inputs['C'],sep='\t')
rho=float(r.loc[r.measurement.eq('DNMT3A_mean'),'Spearman'].iloc[0]);assert len(g)==20 and g.LADs.sum()==1405
assert m.windows.sum()==457584 and len(m)==92 and m.windows.iloc[-1]==2584
assert m.x.is_monotonic_increasing
for font in ['arial.ttf','arialbd.ttf']:font_manager.fontManager.addfont(str(args.fonts/font))
plt.rcParams.update({'font.family':'Arial','font.size':8.5,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.65,'xtick.major.width':.65,'ytick.major.width':.65})
fig=plt.figure(figsize=(170/25.4,205/25.4))
def old_box(left,bottom,width,height): return [left,(bottom*145+60)/205,width,height*145/205]
axA=fig.add_axes(old_box(.095,.625,.365,.285));axB=fig.add_axes(old_box(.615,.625,.365,.285));axC=fig.add_axes(old_box(.095,.105,.885,.34))
axA.plot(g.B1,g.DNMT3A_mean,c='#416f98',lw=1.5)
axA.set(xlim=(-.04,1.94),ylim=(-.51,.17),xticks=[0,.5,1,1.5],yticks=[-.4,-.2,0],xlabel=r'Culture TTX log$_2$(B1/Input) per LAD',ylabel=r'Mean log$_2$(DNMT3A/Input)')
axA.set_title('DNMT3A in postnatal cortex, 2 weeks\n(Stroud et al., 2017)',fontsize=9,pad=9)
axA.text(.05,.09,f'Spearman ρ = {rho:.2f}',transform=axA.transAxes,fontsize=8.5)
axB.axvspan(0,1,color='#eeeeee',zorder=0)
for x in [0,1]:axB.axvline(x,c='#a0a0a0',lw=.55,zorder=1)
axB.plot(b.x,b.H3K9me2_H3mean,c='#a96e28',lw=1.5,zorder=2)
axB.set(xlim=(-1.12,2.12),ylim=(-.82,.25),xticks=[-1,0,.5,1,2],xticklabels=['−500 kb','Start','LAD','End','+500 kb'],yticks=[-.6,-.3,0],xlabel='Culture TTX LAD body and flanks',ylabel=r'Mean log$_2$(H3K9me2/H3)')
axB.set_title('H3K9me2 in postnatal cortex, 2 weeks\n(Stroud et al., 2017)',fontsize=9,pad=9)
for column,label,color in [('1wk','Postnatal cortex, 1 week (Stroud et al., 2017)','#292929'),('2wk','Postnatal cortex, 2 weeks (Stroud et al., 2017)','#3f51b5'),('Baseline','Culture baseline','#25836b')]:
 axC.plot(m.x,m[column],lw=1.5,color=color,label=label)
axC.set(xlim=(-1.75,1.30),xticks=[-1.5,-1,-.5,0,.5,1],ylim=(45,90),yticks=[50,60,70,80,90],xlabel=r'Postnatal cortex, 2 weeks — mean log$_2$(DNMT3A/Input)',ylabel='mCG (%)')
axC.set_title('mCG across DNMT3A enrichment',fontsize=9,pad=9)
axC.legend(frameon=False,loc='lower right',fontsize=8.5,handlelength=2.5,labelspacing=.35)
for letter,x,y in [('A',.018,.945),('B',.533,.945),('C',.018,.480)]:fig.text(x,(y*145+60)/205,letter,fontweight='bold',fontsize=12,va='center')
trends=pd.read_csv(inputs['DE_trends'],sep='\t',float_precision='round_trip')
stats=pd.read_csv(inputs['DE_correlations'],sep='\t',float_precision='round_trip')
axD=fig.add_axes([.095,14.35/205,.365,40.5/205]);axE=fig.add_axes([.615,14.35/205,.365,40.5/205])
de_metrics={}
for ax,scope,title,n in [(axD,'All eligible bins','Genome-wide B1–mCG',78007),(axE,'TTX LAD bins','B1–mCG within LADs',35676)]:
    q=trends.loc[trends.scope.eq(scope)].copy()
    result=stats.loc[stats.scope.eq(scope)].iloc[0]
    assert len(q)==20 and int(q.n.sum())==int(result.n_bins)==n
    x=q.x.to_numpy();y=q['median'].to_numpy();lo=q.lo.to_numpy();hi=q.hi.to_numpy()
    ax.fill_between(x,lo,hi,color='#88b7c1',alpha=.35,linewidth=.6)
    ax.plot(x,y,color='#356c77',lw=1.5,marker='o',markersize=2.5,markeredgewidth=0)
    ax.set(xlim=(-1,2),ylim=(55,95),xticks=[-1,0,1,2],yticks=[55,65,75,85,95],xlabel=r'Culture TTX log$_2$(B1/Input)',ylabel='Culture baseline mCG (%)')
    ax.set_title(title,fontsize=9,pad=9)
    value=float(result.Spearman_rho)
    label=f'Spearman ρ = {value:+.2f}'.replace('-', '−')
    ax.text(.05,.9,label,transform=ax.transAxes,fontsize=8.5)
    np.testing.assert_array_equal(ax.lines[0].get_ydata(),y)
    assert np.min(lo)>=55 and np.max(hi)<=95
    de_metrics[scope]={'n_windows':n,'Spearman_rho':value,'groups':len(q)}
for letter,x in [('D',.018),('E',.533)]:fig.text(x,60/205,letter,fontweight='bold',fontsize=12,va='center')

# Check plotted coordinates and axis limits.
for ax in [axA,axB,axC,axD,axE]:
 for line in ax.lines:
  if len(line.get_xdata())<=2:continue
  xx=np.asarray(line.get_xdata(),float);yy=np.asarray(line.get_ydata(),float)
  assert xx.min()>=ax.get_xlim()[0] and xx.max()<=ax.get_xlim()[1]
  assert yy.min()>=ax.get_ylim()[0] and yy.max()<=ax.get_ylim()[1]
np.testing.assert_array_equal(axA.lines[0].get_ydata(),g.DNMT3A_mean)
np.testing.assert_array_equal(axB.lines[-1].get_ydata(),b.H3K9me2_H3mean)
for line,col in zip(axC.lines,['1wk','2wk','Baseline']):np.testing.assert_array_equal(line.get_ydata(),m[col])
fig.canvas.draw();renderer=fig.canvas.get_renderer();fw,fh=fig.canvas.get_width_height()
for t in fig.findobj(matplotlib.text.Text):
 if t.get_visible() and t.get_text():
  box=t.get_window_extent(renderer);assert box.x0>=-1 and box.y0>=-1 and box.x1<=fw+1 and box.y1<=fh+1,(t.get_text(),box)
name='Fig_S10';outputs=[]
for ext in ['png','pdf','svg']:
 dst=FO/(name+'.'+ext);fig.savefig(dst,dpi=300);outputs.append(dst)
plt.close(fig)
