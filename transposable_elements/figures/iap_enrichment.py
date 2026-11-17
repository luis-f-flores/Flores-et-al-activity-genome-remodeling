"""IAP hypo-DMR enrichment and genomic-context figure panels."""
import csv
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pathlib import Path
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'iap_matching')
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--fonts',type=Path,default=None)
args=ap.parse_args(); IN=args.inputs; OUT=args.output; OUT.mkdir(parents=True,exist_ok=True)

font_dir=args.fonts
if font_dir is None:
    font_dir=next((p for p in [Path('C:/Windows/Fonts'),Path('/mnt/c/Windows/Fonts'),Path('/usr/share/fonts/truetype/msttcorefonts')] if (p/'arial.ttf').exists()),None)
if font_dir is None: raise SystemExit('Set --fonts to the directory containing Arial .ttf files.')

pdfmetrics.registerFont(TTFont('Arial',str(font_dir/'arial.ttf')))
pdfmetrics.registerFont(TTFont('ArialBold',str(font_dir/'arialbd.ttf')))
rows=list(csv.DictReader((IN/'genomic_enrichment_results.tsv').open(),delimiter='	'))
assert [(int(x['n']),int(x['observed'])) for x in rows]==[(61,9),(38,6)]
ORANGE='#d95f02';GRAY='#b6b6b6';PURPLE='#7b3294'
def text(c,x,y,s,size=7.5,font='Arial',align='left',color='black'):
    c.setFont(font,size);c.setFillColor(color)
    getattr(c,{'left':'drawString','center':'drawCentredString','right':'drawRightString'}[align])(x,y,s)
def line(c,x1,y1,x2,y2,color='black',width=.6):
    c.setStrokeColor(color);c.setLineWidth(width);c.line(x1,y1,x2,y2)
def rect(c,x,y,w,h,col):c.setFillColor(col);c.rect(x,y,w,h,stroke=0,fill=1)

# Figure 7F: observed hypo-DMR overlap frequency and matched genomic expectation.
c=canvas.Canvas(str(OUT/'Fig_7F.pdf'),pagesize=(122,185))
text(c,31,172,'F',11,'ArialBold')
x0,y0,pw,ph=31,42,86,106;ymax=22
for tick in (0,5,10,15,20):
    yy=y0+ph*tick/ymax;line(c,x0,yy,x0+pw,yy,'#dddddd',.4);text(c,x0-4,yy-2.7,str(tick),7.5,align='right')
for center,row,label in zip([51,96],rows,['Fragments','Grouped']):
    obs=float(row['observed_pct']);exp=float(row['expected_pct']);bw=13
    rect(c,center-bw,y0,bw,ph*obs/ymax,ORANGE);rect(c,center,y0,bw,ph*exp/ymax,GRAY)
    text(c,center-6.5,y0+ph*obs/ymax+11,f'{obs:.1f}%',7,align='center')
    text(c,center-6.5,y0+ph*obs/ymax+3,f"{row['observed']}/{row['n']}",7,align='center')
    text(c,center+6.5,y0+ph*exp/ymax+3,f'{exp:.1f}%',7,align='center')
    text(c,center,157,f"{float(row['fold']):.2f}-fold",7.3,align='center')
    text(c,center,148,f"P = {float(row['p_enrichment']):.4f}",7,align='center')
    text(c,center,y0-12,label,7.3,align='center')
line(c,x0,y0,x0,y0+ph);line(c,x0,y0,x0+pw,y0)
c.saveState();c.translate(9,96);c.rotate(90);text(c,0,0,'Hypo-DMR overlap (%)',8.2,align='center');c.restoreState()
rect(c,31,14,8,5,ORANGE);text(c,43,13,'Observed',7.1)
rect(c,31,4,8,5,GRAY);text(c,43,3,'Matched genomic',7.1)
c.save()

# Set the S17F panel height while retaining its text size.
c=canvas.Canvas(str(OUT/'Fig_S17F.pdf'),pagesize=(128,124))
text(c,10,112,'F',11,'ArialBold');text(c,76,112,'Gene-body context',7.4,align='center')
text(c,76,101,'Fisher P = 0.725',7.1,align='center')
x0,y0,pw,ph=28,25,91,69
for tick in (0,10,20):
    yy=y0+ph*tick/25;line(c,x0,yy,x0+pw,yy,'#dddddd',.4);text(c,x0-4,yy-2.5,str(tick),7,align='right')
for xx,num,den,col,lab in [(51,6,36,PURPLE,'Inside'),(96,3,25,GRAY,'Outside')]:
    v=100*num/den;rect(c,xx-12,y0,24,ph*v/25,col);text(c,xx,y0+ph*v/25+10,f'{v:.1f}%',7,align='center');text(c,xx,y0+ph*v/25+2,f'{num}/{den}',7,align='center');text(c,xx,14,lab,7.2,align='center')
text(c,74,3,'gene bodies',7.2,align='center');line(c,x0,y0,x0,y0+ph);line(c,x0,y0,x0+pw,y0)
c.saveState();c.translate(8,63);c.rotate(90);text(c,0,0,'Hypo-DMR overlap (%)',7.4,align='center');c.restoreState();c.save()

