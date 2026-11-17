"""Draw gene-body methylation–RNA correlations and fitted regression lines."""
from pathlib import Path
from io import BytesIO
import argparse, json, math, hashlib
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
import pypdfium2 as P

ap=argparse.ArgumentParser(description="Draw gene-body methylation–RNA relationships.")
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'gene_body_methylation')
ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'output')
ap.add_argument('--fonts',type=Path,default=Path('C:/Windows/Fonts') if Path('C:/Windows/Fonts').exists() else Path('/mnt/c/Windows/Fonts'))
args=ap.parse_args();args.output.mkdir(parents=True,exist_ok=True)
D=json.loads((args.inputs/'panel_data.json').read_text())
for name,f in [('Arial','arial.ttf'),('ArialBold','arialbd.ttf')]:pdfmetrics.registerFont(TTFont(name,str(args.fonts/f)))
BLACK='#171717';GRAY='#aaaaaa';BLUE='#269dcc';RED='#f14572'
def text(c,x,y,s,size=8.5,align='left',color=BLACK,rotate=0):
    c.saveState();c.setFillColor(HexColor(color));c.setFont('Arial',size);c.translate(x,y);c.rotate(rotate)
    getattr(c,{'left':'drawString','center':'drawCentredString','right':'drawRightString'}[align])(0,0,s);c.restoreState()

def line(c,xy,color=BLACK,width=.5,dash=None):
    c.saveState();c.setStrokeColor(HexColor(color));c.setLineWidth(width)
    if dash:c.setDash(*dash)
    p=c.beginPath();p.moveTo(*xy[0]);[p.lineTo(*a) for a in xy[1:]];c.drawPath(p);c.restoreState()

def polygon(c,xy,color,alpha=1,stroke=None,width=.5):
    c.saveState();c.setFillColor(HexColor(color));c.setFillAlpha(alpha)
    if stroke:c.setStrokeColor(HexColor(stroke));c.setLineWidth(width)
    p=c.beginPath();p.moveTo(*xy[0]);[p.lineTo(*a) for a in xy[1:]];p.close();c.drawPath(p,fill=1,stroke=bool(stroke));c.restoreState()

class Ax:
    def __init__(self,c,rect,xlim,ylim):self.c=c;self.x0,self.y0,self.w,self.h=rect;self.xlim=xlim;self.ylim=ylim
    def xy(self,x,y):return (self.x0+(x-self.xlim[0])/(self.xlim[1]-self.xlim[0])*self.w,self.y0+(y-self.ylim[0])/(self.ylim[1]-self.ylim[0])*self.h)
    def path(self,x,y,color=BLACK,width=.7,dash=None):line(self.c,[self.xy(a,b) for a,b in zip(x,y)],color,width,dash)
    def axes(self,xt,yt,ylabels=True,xlabels=True,xlab=None,ylab=None,ysize=8.5):
        line(self.c,[(self.x0,self.y0+self.h),(self.x0,self.y0),(self.x0+self.w,self.y0)])
        for v in xt:
            x,y=self.xy(v,self.ylim[0]);line(self.c,[(x,y),(x,y-2)])
            if xlabels:text(self.c,x,y-10,f'{v:g}'.replace('-','−'),8,align='center')
        for v in yt:
            x,y=self.xy(self.xlim[0],v);line(self.c,[(x,y),(x-2,y)])
            if ylabels:text(self.c,x-4,y-2.5,f'{v:g}'.replace('-','−'),8,align='right')
        if xlab:text(self.c,self.x0+self.w/2,self.y0-22,xlab,9,align='center')
        if ylab:text(self.c,self.x0-20,self.y0+self.h/2,ylab,ysize,align='center',rotate=90)
    def clip(self):
        self.c.saveState();p=self.c.beginPath();p.rect(self.x0,self.y0,self.w,self.h);self.c.clipPath(p,stroke=0)
    def unclip(self):self.c.restoreState()
    def zero(self,vertical=True,horizontal=True,dash=None):
        if horizontal and self.ylim[0]<0<self.ylim[1]:self.path(self.xlim,[0,0],GRAY,.5,dash)
        if vertical and self.xlim[0]<0<self.xlim[1]:self.path([0,0],self.ylim,GRAY,.5,dash)
    def scatter(self,x,y,color=BLUE,r=.95,alpha=.8):
        assert all(self.xlim[0]<=v<=self.xlim[1] for v in x),'Point outside x axis'
        assert all(self.ylim[0]<=v<=self.ylim[1] for v in y),'Point outside y axis'
        self.c.saveState();self.c.setFillColor(HexColor(color));self.c.setFillAlpha(alpha)
        for xx,yy in zip(x,y):self.c.circle(*self.xy(xx,yy),r,stroke=0,fill=1)
        self.c.restoreState()
    def band(self,x,lo,hi,color,alpha=.15):
        polygon(self.c,[self.xy(a,b) for a,b in zip(x,lo)]+[self.xy(a,b) for a,b in zip(x[::-1],hi[::-1])],color,alpha)

def pformat(p):return f'{p:.3f}' if p>=.001 else f'{p:.3g}'

def draw_gene_body_panels(c):
    for k,cohort in enumerate(['differentially_methylated','DEG_subset']):
        for j,tp in enumerate([1,6,24]):
            d=next(x for x in D['S06'] if x['cohort']==cohort and x['tp']==tp);x=np.array(d['x']);y=np.array(d['y'])
            pad=(max(x)-min(x))*.05;xl=(min(x)-pad,max(x)+pad)
            # Keep all measurements in view, with a shared symmetric RNA axis in each row.
            row=[z for z in D['S06'] if z['cohort']==cohort];ym=max(3,math.ceil(max(abs(v) for z in row for v in z['y'])+.1))
            a=Ax(c,(51.0+161.6*j,169 if k==0 else 33,122,90.7),xl,(-ym,ym))
            xt=[v for v in ([-25,0,25] if k==0 else ([0,5,10] if tp==1 else [-10,0,10])) if xl[0]<=v<=xl[1]]
            a.axes(xt,[-ym,0,ym],xlab='ΔmCG (%)' if k else None,ylab='RNA log2 change' if j==0 else None)
            a.clip();a.scatter(x,y,'#518bc8',.9,.34);a.zero(vertical=False);f=d['fit'];a.band(f['x'],f['lo'],f['hi'],'#e2564f',.16);a.path(f['x'],f['y'],'#ed514b',.8);a.unclip()
            text(c,a.x0+3,a.y0+a.h-8,f"r = {d['r']:.3f}; P = {d['p']:.3g}",8);text(c,a.x0+3,a.y0+a.h-17,f"n = {d['n']:,}",8)
            if not k:text(c,a.x0+a.w/2,265,['1 h','6 h','Recovery'][j],9,align='center')
    text(c,51.6,142,'DEG subset',8.5)

def render(p,dpi):
    d=P.PdfDocument(p);pg=d[0];bm=pg.render(scale=dpi/72);im=bm.to_pil().convert('RGB').copy();bm.close();pg.close();d.close();return im

# Draw the quantitative panels on their own canvas.
pdf=args.output/'Fig_S14.pdf'
c=canvas.Canvas(str(pdf),pagesize=(530,290),invariant=1)
draw_gene_body_panels(c)
c.save()
render(pdf,300).save(pdf.with_suffix('.png'))
print(pdf)
