"""Draw LAD-segment B1–RNA correlations and RNA changes by B1 category."""
from pathlib import Path
from io import BytesIO
import argparse, json, math, hashlib
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
import pypdfium2 as P

ap=argparse.ArgumentParser(description="Draw B1–RNA relationships.")
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'b1_rna')
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

def box(c,a,x,d,col):
    q1,med,q3=d['q25'],d['median'],d['q75'];xy=[a.xy(x-.27,q1),a.xy(x+.27,q1),a.xy(x+.27,q3),a.xy(x-.27,q3)];polygon(c,xy,col,1,BLACK,.6)
    a.path([x-.27,x+.27],[med,med],BLACK,.8)
    for lo,hi in [(d['whislo'],q1),(q3,d['whishi'])]:a.path([x,x],[lo,hi],BLACK,.6)
    for y in [d['whislo'],d['whishi']]:a.path([x-.14,x+.14],[y,y],BLACK,.7)

def draw_b1_rna_panels(c):
    for j,d in enumerate(D['S17']):
        a=Ax(c,(45.35+170.08*j,305.3,113.4,122),(-1.3,1.3),(-7,7));a.axes([-1,0,1],[-6,-4,-2,0,2,4,6],ylabels=j==0)
        if j==0:text(c,16,a.y0+a.h/2,'RNA log2 change',9,align='center',rotate=90);text(c,27,a.y0+a.h/2,'(region median)',8.5,align='center',rotate=90)
        a.clip();a.zero(dash=(1.5,1.5));a.scatter(d['x'],d['y'],'#71868e',.85,.4);a.unclip();text(c,a.x0+3,a.y0+a.h-9,f"ρ = {d['rho']:.2f}",8);text(c,a.x0+3,a.y0+a.h-18,f"P = {d['p']:.4f}",8)
    text(c,266.3,278,'ΔB1 (log2 B1/Input)',9,align='center')
    for j,d in enumerate(D['S17boxes']):
        a=Ax(c,(59.5+255.1*j,65.3,181.4,101.5),(.2,3.8),(-3,2.5));a.axes([],[-3,-2,-1,0,1,2],ylab='RNA log2 change' if j==0 else None)
        a.clip();a.zero(vertical=False,dash=(1.5,1.5))
        for i,(v,col) in enumerate(zip(d['groups'],['#bb72a6','#cbd0d3','#008775'])):box(c,a,i+1,v,col)
        a.path([1,1,3,3],[2.15,2.3,2.3,2.15],'#555555',.5);a.unclip();text(c,a.xy(2,0)[0],a.xy(2,2.37)[1],f"P = {d['p']:.3g}",8,align='center')
        for i,v in enumerate(d['groups']):
            x=a.xy(i+1,0)[0];lab=['Lower B1','Change not','Higher B1'][i];text(c,x,54,lab,8,align='center')
            if i==1:text(c,x,46,'detected',8,align='center');ny=37
            else:ny=45
            text(c,x,ny,f"n = {v['n']:,}",8,align='center')

def render(p,dpi):
    d=P.PdfDocument(p);pg=d[0];bm=pg.render(scale=dpi/72);im=bm.to_pil().convert('RGB').copy();bm.close();pg.close();d.close();return im

# Draw the quantitative panels on their own canvas.
pdf=args.output/'Fig_S09.pdf'
c=canvas.Canvas(str(pdf),pagesize=(530,455),invariant=1)
draw_b1_rna_panels(c)
c.save()
render(pdf,300).save(pdf.with_suffix('.png'))
print(pdf)
