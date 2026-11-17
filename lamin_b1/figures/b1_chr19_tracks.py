"""Draw the chromosome 19 B1 tracks for Figure 1C."""
import os
os.environ.update(MPLBACKEND='Agg',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
from pathlib import Path
import argparse
import numpy as np
P=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'b1_chr19'
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--inputs',type=Path,default=P);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--font-dir',type=Path,required=True)
args=ap.parse_args();P=args.inputs;args.assemble_only=False
O=args.out.resolve();O.mkdir(parents=True,exist_ok=True)
genes=[line.split('\t') for line in (P/'genes.bed').read_text().splitlines()]
data=dict(np.load(P/'chr19_tracks.npz'))
starts,sig,means,state,ws,wm,length=[data[k] for k in ['starts','sig','means','state','ws','wm','length']]
length=int(length);oy=(float(wm.min()),float(wm.max()))
np.testing.assert_array_equal(starts,np.arange(10800000,21600000,30000))
IPS=['L9','L12','L11','L14','L16','L10','L13','L15'];GROUPS=[[0,1],[2,3,4],[5,6,7]]
COL=['#0072B2','#D55E00','#4D4D4D'];SC=['#e1e4e8','#78b9db','#154968']
component=O/'C_tracks.pdf'
if not args.assemble_only:
    from matplotlib import pyplot as plt,font_manager
    from matplotlib.patches import Rectangle
    for n in ['arial.ttf','arialbd.ttf']:font_manager.fontManager.addfont(str(args.font_dir/n))
    plt.rcParams.update({'font.family':'Arial','font.size':8.5,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'xtick.labelsize':8.5,'ytick.labelsize':8.5})
    PW,PH=510.236,590.161
    fig=plt.figure(figsize=(PW/72,PH/72));fig.patch.set_alpha(0)
    def ax(x,y,w,h):return fig.add_axes([x/PW,y/PH,w/PW,h/PH])
    def txt(x,y,t,size=8.5,**kw):
        kw.setdefault('va','baseline')
        return fig.text(x/PW,y/PH,t,fontsize=size,**kw)
    def segments(ix,ss):return np.split(ix,np.flatnonzero(np.diff(ss[ix])!=30000)+1)
    curvechecks=[]
    for j,(g,x0) in enumerate(zip(GROUPS,[113.3858,240.9449,368.5039])):
        width=110.5513
        aa=ax(x0,462.348,width,8.086)
        for ix in segments(np.arange(len(ws)),ws):aa.plot((ws[ix]+15000)/1e6,wm[ix,j],color=COL[j],lw=.45)
        aa.add_patch(Rectangle((10.8,oy[0]),10.8,oy[1]-oy[0],facecolor='#b6a7cf',edgecolor='#8e86a2',alpha=.3,lw=.4,zorder=-1))
        aa.set(xlim=(0,length/1e6),ylim=oy);aa.axis('off')
        aa=ax(x0,413.2414,width,45.354)
        # Draw the condition mean as a signed area.
        x=(starts+15000)/1e6
        y=means[:,j]
        aa.fill_between(x,0,y,where=y>=0,color=COL[j],interpolate=True,linewidth=0)
        aa.fill_between(x,0,y,where=y<=0,color='#999999',interpolate=True,linewidth=0)
        np.testing.assert_array_equal(y,data['means'][:,j]);curvechecks.append(True)
        aa.axhline(0,color='#aaa',lw=.35,zorder=-2)
        aa.set(xlim=(10.8,21.6),ylim=(-1,2.5),yticks=[-1,0,1,2],xticks=[])
        aa.spines['bottom'].set_visible(False);aa.tick_params(length=2,pad=1)
        aa=ax(x0,401.903,width,5.669)
        for s,k in zip(starts,state[:,j]):aa.add_patch(Rectangle((s/1e6,0),.03,1,facecolor=SC[k],lw=0))
        aa.set(xlim=(10.8,21.6),ylim=(0,1));aa.axis('off')
        for dr,y in [('gained',392.719),('lost',385.575)]:
            aa=ax(x0,y,width,3.968)
            if j:
                mask=(state[:,j]>0)&(state[:,0]==0) if dr=='gained' else (state[:,j]==0)&(state[:,0]>0)
                for s in starts[mask]:aa.add_patch(Rectangle((s/1e6,0),.03,1,facecolor='#009e73' if dr=='gained' else '#cc79a7',lw=0))
            aa.set(xlim=(10.8,21.6),ylim=(0,1));aa.axis('off')
        # Compact, coordinate-faithful gene spans, separated by strand.
        aa=ax(x0,375.8,width,5.2)
        for chrom,gs,ge,name,score,strand in genes:
            l=max(int(gs),10800000)/1e6;r=min(int(ge),21600000)/1e6
            yy=.72 if strand=='+' else .24
            aa.plot([l,r],[yy,yy],color='#263e54',lw=.65,solid_capstyle='butt')
            if r-l>.08:
                mid=(l+r)/2;xx=mid+.017 if strand=='+' else mid-.017
                aa.plot([xx-.017 if strand=='+' else xx+.017,xx,xx-.017 if strand=='+' else xx+.017],[yy-.12,yy,yy+.12],color='#263e54',lw=.4)
        aa.set(xlim=(10.8,21.6),ylim=(0,1));aa.axis('off')
        aa=ax(x0,371.2257,width,.1);aa.set(xlim=(10.8,21.6),ylim=(0,1),yticks=[],xticks=[10.8,16.2,21.6],xticklabels=['10.8','16.2','21.6'])
        aa.spines[['left','top','right']].set_visible(False);aa.tick_params(length=2,pad=1)
        aa.get_xticklabels()[0].set_ha('left');aa.get_xticklabels()[-1].set_ha('right')
    txt(21.6702,466.35,'chr19')
    txt(21.6702,457.25,'0-61 Mb',7)
    txt(21.5937,364.65,'chr19 (Mb)')
    txt(21.5937,375.4,'Genes')
    txt(88.5,435.92,'log2(B1/Input)',rotation=90,ha='center',va='center')
    txt((113.3858+223.9371)/2,389.1,'Reference',ha='center',color='#777')
    component=O/'C_tracks.pdf';fig.savefig(component,transparent=True,metadata={'CreationDate':None,'ModDate':None});plt.close(fig)


print('Created '+str(component))
