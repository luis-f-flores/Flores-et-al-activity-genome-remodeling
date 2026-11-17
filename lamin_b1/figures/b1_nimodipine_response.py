"""Draw B1 and RNA responses to nimodipine for the 80 LAD-associated DEGs."""
from pathlib import Path
import argparse,csv
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--fonts',type=Path,required=True)
args=p.parse_args();args.output.mkdir(parents=True,exist_ok=False)
details=args.inputs/'fixed_Bic_DEG_return_details_125.tsv'
with details.open() as f:
    rows = [r for r in csv.DictReader(f, delimiter="\t") if r["LAD_Bic_vs_TTX"] == "True"]
assert len(rows) == len({r["gene_id"] for r in rows})
counts = [0, 0, 0, 0]
for r in rows:
    b, bn, n = (float(r[k]) for k in ["b", "bn", "n"])
    f, fn, direct = (float(r[k]) for k in ["fc_Bic_vs_TTX", "fc_Bic_Nimo_vs_TTX", "fc_Bic_Nimo_vs_Bic"])
    bt = b*n < 0 and abs(bn) < abs(b)
    rt = f*direct < 0 and abs(fn) < abs(f)
    assert bt == (r["B1_toward"] == "True")
    assert rt == (r["RNA_toward"] == "True")
    counts[0 if bt and rt else 1 if bt else 2 if rt else 3] += 1


for name,fn in [('Arial','arial.ttf'),('ArialBold','arialbd.ttf')]:pdfmetrics.registerFont(TTFont(name,str(args.fonts/fn)))
mm=72/25.4;extra=66*mm
c=canvas.Canvas(str(args.output/'Fig_S09F.pdf'),pagesize=(170*mm,extra))
c.setFillColor(HexColor("#171717"))
c.setFont("ArialBold", 9)
c.drawString(2*mm, extra-6*mm, "F")
c.setFont("Arial", 9)
c.drawString(10*mm, extra-6*mm, "Movement toward TTX after nimodipine")
c.setFont("Arial", 8)
c.setFillColor(HexColor("#555555"))
c.drawString(10*mm, extra-12*mm, f"{len(rows)} bicuculline-responsive LAD genes")
x0, x1, axis_y = 53*mm, 155*mm, 12*mm
scale = (x1-x0)/40
ys = [(extra/mm-23-i*8)*mm for i in range(4)]
colors = ["#65539D", "#3F8F9D", "#C3834F", "#B7BDC5"]
labels = ["Both B1 and RNA", "B1 only", "RNA only", "Neither"]
for label, n, y, color in zip(labels, counts, ys, colors):
    c.setFont("Arial", 8)
    c.setFillColor(HexColor("#171717"))
    c.drawRightString(x0-3*mm, y-2.5, label)
    c.setFillColor(HexColor(color))
    c.rect(x0, y-2*mm, n*scale, 4*mm, fill=1, stroke=0)
    c.setFillColor(HexColor("#171717"))
    c.drawString(x0+n*scale+2*mm, y-2.5, str(n))
c.setStrokeColor(HexColor("#777777"))
c.setLineWidth(.45)
c.line(x0, axis_y, x1, axis_y)
for n in range(0, 41, 10):
    x = x0+n*scale
    c.line(x, axis_y, x, axis_y-1*mm)
    c.setFont("Arial", 7.5)
    c.drawCentredString(x, axis_y-4.5*mm, str(n))
c.setFont("Arial", 8)
c.drawCentredString((x0+x1)/2, 2.5*mm, "Genes")

c.save()
print(f'{len(rows)} genes:',counts)
