'Plot WGBS library coverage and RNA-cluster GO results.\n\nCoverage input is a tab-delimited table with sample, min, max, mean, median,\npct5, pct10 and pct20 columns. GO inputs are five CSV exports named\nCluster_<1..5>_GO_Biological_Process_Top15.xlsx.csv.'
from pathlib import Path
import argparse
import csv
import math
import textwrap
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors

MM = 72 / 25.4

def txt(c,x,y,s,size=9,bold=False,align='left',color='#222222',italic=False):
    c.setFont('Arial-Italic' if italic else 'Arial-Bold' if bold else 'Arial',size)
    c.setFillColor(colors.HexColor(color))
    fn={'left':c.drawString,'center':c.drawCentredString,'right':c.drawRightString}[align]
    fn(x*MM,c._pagesize[1]-y*MM-size*.8,s)

def panel(c,letter,title,x,y):
    txt(c,x,y,letter,11,True)
    txt(c,x+6,y+.4,title,9)

def start(path,w,h):
    c=canvas.Canvas(str(path),pagesize=(w*MM,h*MM),pageCompression=1)
    c.setTitle(path.stem.replace('_',' '))
    return path,c

def write_table(path,rows):
    with path.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

def coverage_plot(source,output):
    with source.open(encoding='utf-8-sig',newline='') as handle:
        rows=list(csv.DictReader(handle,delimiter='\t'))
    assert len(rows)==21 and len({r['sample'] for r in rows})==21
    order=['0','1M','1K','6M','6K','24M','24K']
    labels=['0 h','1 h M','1 h K','6 h M','6 h K','24 h M','24 h K']
    for r in rows:
        for k in ['min','max','mean','median','pct5','pct10','pct20']:
            r[k]=float(r[k])
    out,c=start(output/'Fig_S11A.pdf',180,105)
    columns=['sample','min','max','mean','median','pct5','pct10','pct20']
    labels=['Sample','Min','Max','Mean','Median','≥5× (%)','≥10× (%)','≥20× (%)']
    xs=[4,29,49,69,92,114,137,160]
    for x,label in zip(xs,labels):txt(c,x,4,label,8,True)
    for i,row in enumerate(rows):
        for x,key in zip(xs,columns):
            value=row[key]
            label=str(value) if key=='sample' else f'{value:.2f}' if key in ['mean','pct5','pct10','pct20'] else f'{value:.0f}'
            txt(c,x,10+i*4.2,label,8)
    c.save()
    write_table(output/'wgbs_coverage.tsv',rows)

def go_plot(go_root,output):
    # Panels B–F share legends; the upper-left position is reserved for panel A.
    go=[]
    for cl in range(1,6):
     p=go_root/f'Cluster_{cl}_GO_Biological_Process_Top15.xlsx.csv'
     with p.open(encoding='utf-8-sig',newline='') as handle:
      rr=sorted(list(csv.DictReader(handle)),key=lambda r:float(r['p.adjust']))[:5]
     assert len(rr)==5
     for r in rr:r['cluster']=cl;r['ratio']=int(r['GeneRatio'].split('/')[0])/int(r['GeneRatio'].split('/')[1]);r['logq']=-math.log10(float(r['p.adjust']))
     go+=rr
    maxrat=max(r['ratio'] for r in go)*1.1;maxq=max(r['logq'] for r in go);maxcount=max(int(r['Count']) for r in go)
    width,height=170,155
    out,c=start(output/'Fig_S13_BF.pdf',width,height)
    def label(x,y,value,size=9,bold=False,align='left'):
     txt(c,x,y,value,size,bold,align,color='#000000')
    def wrap(value,available):
     lines=[];line=''
     for word in value.split():
      candidate=(line+' '+word).strip()
      if line and pdfmetrics.stringWidth(candidate,'Arial',7.8)>available*MM:
       lines.append(line);line=word
      else:line=candidate
     if line:lines.append(line)
     return lines
    def dotplot(cl,x,y,letter):
     label(x,y,letter,11,True);label(x+6,y+.4,f'RNA cluster {cl}',9)
     ax=x+57;pw=20;centers=[y+8+k*4.8 for k in range(5)];top=y+5.6;bottom=y+30
     c.saveState();c.setStrokeColor(colors.HexColor('#E2E2E2'));c.setLineWidth(.4)
     for yy in centers:c.line(ax*MM,(height-yy)*MM,(ax+pw)*MM,(height-yy)*MM)
     for fraction in [.25,.5,.75]:
      xx=ax+pw*fraction;c.line(xx*MM,(height-top)*MM,xx*MM,(height-bottom)*MM)
     c.setStrokeColor(colors.HexColor('#B7B7B7'));c.setLineWidth(.5)
     c.rect(ax*MM,(height-bottom)*MM,pw*MM,(bottom-top)*MM,stroke=1,fill=0)
     c.restoreState();c.setLineWidth(.4)
     for yy,r in zip(centers,[r for r in go if r['cluster']==cl]):
      lines=wrap(r['Description'],54);assert len(lines)<=2,(r['Description'],lines)
      for j,line in enumerate(lines):label(ax-2,yy-1.6+(j-(len(lines)-1)/2)*2.9,line,7.8,align='right')
      v=min(1,r['logq']/maxq);c.setFillColor(colors.Color(.88-.70*v,.93-.63*v,.97-.42*v));c.setStrokeColor(colors.HexColor('#426385'))
      xx=ax+pw*r['ratio']/maxrat;radius=.65+math.sqrt(int(r['Count'])/maxcount)
      c.circle(xx*MM,(height-yy)*MM,radius*MM,stroke=1,fill=1)
     for value in [0,maxrat/2,maxrat]:label(ax+pw*value/maxrat,bottom+1,f'{value:.2f}',7,align='center')
     label(ax+pw/2,bottom+4.8,'Gene ratio',8,align='center')
    dotplot(1,89,2,'B')
    dotplot(2,2,64,'C');dotplot(3,89,64,'D')
    dotplot(4,2,102,'E');dotplot(5,89,102,'F')
    label(28,145,'−log10(FDR)',8)
    for k in range(50):
     v=k/49;c.setFillColor(colors.Color(.88-.70*v,.93-.63*v,.97-.42*v));c.rect((52+k*.42)*MM,(height-148)*MM,.43*MM,2*MM,stroke=0,fill=1)
    label(52,150,'0',7);label(73,150,f'{maxq:.0f}',7,align='right');label(90,145,'Gene count',8)
    for x,n in [(116,25),(130,75),(144,125)]:
     c.setFillColor(colors.HexColor('#6E91AF'));c.circle(x*MM,(height-146.7)*MM,(.65+math.sqrt(n/maxcount))*MM,stroke=0,fill=1);label(x,150,str(n),7,align='center')
    c.save()
    write_table(output/'rna_cluster_go.tsv',go)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coverage',type=Path)
    parser.add_argument('--go-dir',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--fonts',type=Path,required=True)
    args=parser.parse_args()
    if args.coverage is None and args.go_dir is None:
        parser.error('Supply --coverage, --go-dir, or both.')
    for nm,fn in [('Arial','arial.ttf'),('Arial-Bold','arialbd.ttf'),('Arial-Italic','ariali.ttf')]:
        pdfmetrics.registerFont(TTFont(nm,str(args.fonts/fn)))
    args.output.mkdir(parents=True,exist_ok=False)
    if args.coverage is not None:
        coverage_plot(args.coverage,args.output)
    if args.go_dir is not None:
        go_plot(args.go_dir,args.output)


if __name__=='__main__':
    main()
