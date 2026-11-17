"""Matched DEG and control B1 changes (Figure 3B)."""
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.lines import Line2D
from pathlib import Path
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'matched_b1')
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--fonts',type=Path,default=None)
ap.add_argument('--statistics',type=Path,default=None)
args=ap.parse_args(); IN=args.inputs; OUT=args.output; OUT.mkdir(parents=True,exist_ok=True)

font_dir=args.fonts
if font_dir is None:
    font_dir=next((p for p in [Path('C:/Windows/Fonts'),Path('/mnt/c/Windows/Fonts'),Path('/usr/share/fonts/truetype/msttcorefonts')] if (p/'arial.ttf').exists()),None)
if font_dir is None: raise SystemExit('Set --fonts to the directory containing Arial .ttf files.')

for fn in ['arial.ttf','arialbd.ttf']: fm.fontManager.addfont(str(font_dir/fn))
plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7.8,'ytick.labelsize':7.8,'axes.unicode_minus':False,'pdf.fonttype':42,'ps.fonttype':42,'text.color':'#252525','axes.labelcolor':'#252525'})
C={'down':'#336A9B','up':'#B96332'};GREY='#92989D';INK='#252525'
def density(values,grid,bw):
    return np.exp(-.5*((grid[:,None]-np.asarray(values)[None,:])/bw)**2).mean(1)/(bw*np.sqrt(2*np.pi))
def style(ax):
    for s in ['top','right']:ax.spines[s].set_visible(False)
    for s in ['left','bottom']:ax.spines[s].set_color('#888888')
    ax.tick_params(length=3)

d=pd.read_csv(IN/'paired_differences.tsv',sep='	')
q=pd.read_csv(IN/'reference_statistics.tsv',sep='	').set_index('RNA_direction')
assert d.RNA_direction.value_counts().to_dict()=={'down':56,'up':24}
rng=np.random.default_rng(202611091)
values=np.r_[d.control_B1_delta,d.DEG_B1_delta]
raw_bw=float(values.std(ddof=1)*len(values)**(-.2))
diff_bw=float(d.excess_delta.std(ddof=1)*len(d)**(-.2))
grid=np.linspace(-.65,.65,600)

# A: all paired values, pooled by class. Each control point is a matched-set mean.
fig,ax=plt.subplots(figsize=(245/72,163/72));fig.subplots_adjust(left=.23,right=.985,bottom=.27,top=.77)
fig.text(.01,.943,'B',ha='left',fontsize=11,fontweight='bold')
fig.text(.59,.94,'Bicuculline vs TTX',ha='center',fontsize=9)
ax.axhline(0,color='#666666',lw=1,zorder=0)
densities=[]
for dr in ['down','up']:
    a=d[d.RNA_direction==dr]
    for col in ['control_B1_delta','DEG_B1_delta']:densities.append(density(a[col],grid,raw_bw))
norm=max(v.max() for v in densities)
counter=0
for dr,offset in [('down',0),('up',2.0)]:
    a=d[d.RNA_direction==dr];r=q.loc[dr];x0=offset;x1=offset+.75
    for x,col,color in [(x0,'control_B1_delta',GREY),(x1,'DEG_B1_delta',C[dr])]:
        v=a[col].to_numpy();de=densities[counter];counter+=1
        mask=(grid>=v.min())&(grid<=v.max());yg=grid[mask];width=.22*de[mask]/norm
        ax.fill_betweenx(yg,x-width,x+width,color=color,alpha=.28,lw=0,zorder=2)
        ax.scatter(x+rng.uniform(-.095,.095,len(v)),v,color=color,s=6,alpha=.56,edgecolor='white',linewidth=.2,zorder=3)
        ax.plot([x-.13,x+.13],[v.mean(),v.mean()],color=INK,lw=1.8,zorder=4)
    ax.plot([x0,x0,x1,x1],[.49,.52,.52,.49],color=INK,lw=.9)
    ax.text((x0+x1)/2,.545,f'P = {r.Holm_P_LAD_cluster:.5f}',ha='center',fontsize=7.8)
    name='RNA-down' if dr=='down' else 'RNA-up'
    ax.text((x0+x1)/2,-.43,f'{name} (n = {len(a)})',transform=ax.get_xaxis_transform(),ha='center',fontsize=8,color=C[dr])
ax.set_xlim(-.5,3.25);ax.set_ylim(-.59,.64);ax.set_yticks([-.4,-.2,0,.2,.4,.6])
ax.set_xticks([0,.75,2,2.75],['Matched\ncontrols','DEGs','Matched\ncontrols','DEGs'])
ax.set_ylabel('ΔB1 [log2(B1/Input)]',labelpad=3);style(ax)
fig.canvas.draw()
for txt in fig.findobj(plt.Text):
    if txt.get_visible() and txt.get_text():
        b=txt.get_window_extent(fig.canvas.get_renderer())
        assert b.x0>=-1 and b.y0>=-1 and b.x1<=fig.bbox.width+1 and b.y1<=fig.bbox.height+1,(txt.get_text(),b.bounds)
fig.savefig(OUT/'Fig_3B.pdf');fig.savefig(OUT/'Fig_3B.png',dpi=180);plt.close(fig)
