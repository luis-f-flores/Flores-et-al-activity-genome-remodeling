"""Draw promoter overlap and sampled significance timing from the recorded tables."""
from pathlib import Path
import argparse,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager
from matplotlib.patches import Circle
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--statistics',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--fonts',type=Path,required=True)
args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
o=json.loads((args.inputs/'plot_values.json').read_text());a=o['overlap'];matched=json.loads(args.statistics.read_text())
assert a['shared']==matched['original']['overlap']==84
for n in ['arial.ttf','arialbd.ttf']:fontManager.addfont(str(args.fonts/n))
plt.rcParams.update({'font.family':'Arial','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.65})
width,height=170,70
fig=plt.figure(figsize=(width/25.4,height/25.4))
def ax(x,y,w,h):return fig.add_axes([x/width,1-(y+h)/height,w/width,h/height])
def txt(x,y,s,**kw):kw.setdefault('va','top');return fig.text(x/width,1-y/height,s,**kw)
def panel(x,y,l,title):txt(x,y,l,fontweight='bold',fontsize=11);txt(x+5,y,title,fontsize=9)
panel(2,3,'A','Promoter DMG–DEG overlap');panel(87,3,'B','Relative timing')
va=ax(4,15,79,39);va.set_xlim(0,79);va.set_ylim(0,39);va.set_aspect('equal');va.axis('off')
va.add_patch(Circle((28,19),18,color='#9DB9DE',alpha=.82,lw=0));va.add_patch(Circle((49,19),14,color='#F99B97',alpha=.76,lw=0))
for x,s in [(20,f"{a['DEG_only']:,}"),(39,str(a['shared'])),(54,str(a['DMG_only']))]:va.text(x,19,s,ha='center',va='center',fontsize=9,fontweight='bold' if x==39 else None)
txt(22,13,'All DEGs',ha='center',color='#3366BB',fontsize=8.5);txt(59,13,'All DMGs',ha='center',color='#D74343',fontsize=8.5)
txt(44,61,f"P = {matched['tests'][0]['p_greater']:.3f}",ha='center',fontsize=8)
ta=ax(99,14,65,39);v=[z['n'] for z in o['timing']];assert v==[27,22,35]
ta.bar(range(3),v,color=['#ED746C','#2BAD50','#759BD1'],width=.72)
ta.set_ylim(0,40);ta.set_yticks([0,10,20,30,40]);ta.set_ylabel('Genes',fontsize=9,labelpad=3)
ta.set_xticks(range(3),['Coincident','Expression\nfirst','Methylation\nfirst'],fontsize=7.7);ta.tick_params(axis='both',labelsize=8,length=2,pad=2)
for i,z in enumerate(o['timing']):ta.text(i,z['n']+.8,f"{z['percentage']:.0f}%",ha='center',fontsize=8)
txt(131.5,64,f"Timing P = {a['timing_chisquare_P']:.3f}",ha='center',fontsize=8)

fig.savefig(args.output/'Fig_S15AB.pdf');plt.close(fig)
