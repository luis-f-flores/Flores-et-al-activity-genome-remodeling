"""Draw promoter methylation–RNA correlations, OLS slopes and R-squared summaries."""
from pathlib import Path
from io import BytesIO
import argparse, json, math, hashlib
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
import pypdfium2 as P

ap=argparse.ArgumentParser(description="Draw promoter methylation–RNA relationships.")
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'promoter_methylation')
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

def draw_promoter_panels(c):
    for j,d in enumerate(D['M03']):
        a=Ax(c,(38.6+167.25*j,294.5,126.0,59.7),(-65,50),(-2,4.5));a.axes([-60,-30,0,30],[-2,0,2,4],xlab='ΔmCG (%)',ylab='RNA log2 change' if j==0 else None)
        a.clip();a.zero();f=d['fit'];a.band(f['x'],f['lo'],f['hi'],'#777777',.16);a.path(f['x'],f['y'],'#555555',.7,dash=(2,1))
        x=np.asarray(d['x']);y=np.asarray(d['y'])
        for v,col in [(x<0,BLUE),(x>=0,RED)]:a.scatter(x[v],y[v],col,1.15,.9)
        a.unclip();text(c,a.x0+3,a.y0+a.h-8,f"r = {d['r']:.2f}",8);text(c,a.x0+3,a.y0+a.h-17,f"P = {pformat(d['p'])}",8)
        title=f"{['1 h','6 h','Recovery'][j]} · n = {d['n']}";text(c,a.x0+a.w/2,361,title,9,align='center')
    cols=['#7c7979','#e66c5d','#aaa8a8'];xs=[1,2,3]
    a=Ax(c,(39,177,84,55),(.5,3.5),(-1,1));a.axes([],[-1,-.5,0,.5,1],ylab='OLS slope ± SE');a.zero(vertical=False)
    b=Ax(c,(161,177,84,55),(.5,3.5),(0,.15));b.axes([],[0,.05,.1,.15],ylab='R²')
    for j,d in enumerate(D['M03']):
        x=xs[j];lo,hi=d['slope']-d['se'],d['slope']+d['se'];p=a.xy(x,d['slope']);p0=a.xy(x,0);left=a.xy(x-.36,0)[0];right=a.xy(x+.36,0)[0]
        polygon(c,[(left,p0[1]),(right,p0[1]),(right,p[1]),(left,p[1])],cols[j]);a.path([x,x],[lo,hi],BLACK,.6)
        for y in [lo,hi]:a.path([x-.1,x+.1],[y,y],BLACK,.6)
        py=a.xy(x,max(0,hi))[1]+5 if d['slope']<0 else a.xy(x,hi)[1]+6
        text(c,p[0],py,pformat(d['p']),7.8,align='center')
        y=b.xy(x,d['R2'])[1];base=b.xy(x,0)[1];l=b.xy(x-.42,0)[0];r=b.xy(x+.42,0)[0]
        polygon(c,[(l,base),(r,base),(r,y),(l,y)],'#8272ab');text(c,(l+r)/2,y+4,pformat(d['p']),7.8,align='center')
        for q in [a,b]:text(c,q.xy(x,0)[0],166,['1 h','6 h','Recovery'][j],8,align='center')

def render(p,dpi):
    d=P.PdfDocument(p);pg=d[0];bm=pg.render(scale=dpi/72);im=bm.to_pil().convert('RGB').copy();bm.close();pg.close();d.close();return im

# Draw the quantitative panels on their own canvas.
pdf=args.output/'Fig_6.pdf'
c=canvas.Canvas(str(pdf),pagesize=(530,390),invariant=1)
draw_promoter_panels(c)
c.save()
render(pdf,300).save(pdf.with_suffix('.png'))
print(pdf)
