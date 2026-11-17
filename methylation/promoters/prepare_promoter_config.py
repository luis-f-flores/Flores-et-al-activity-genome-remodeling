"""Write input paths for promoter overlap and conditional matching."""
from pathlib import Path
import argparse,json
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--features',type=Path,required=True);p.add_argument('--master-fpkm',type=Path,required=True)
p.add_argument('--rna',type=Path,nargs=3,required=True);p.add_argument('--windows',type=Path,nargs=3,required=True)
p.add_argument('--promoters',type=Path,nargs=6,required=True,help='Hyper/hypo tables at 1 h, 6 h and recovery.')
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();path=lambda p:str(p.resolve())
c={'package_root':'.','promoter_features':path(a.features),'master_fpkm':path(a.master_fpkm),'deg_padj':.01,'deg_fold_change':1.3,'dmg_qvalue':.01,'dmg_delta':25,'timepoints':{}}
for i,t in enumerate([1,6,24]):c['timepoints'][str(t)]={'rna':path(a.rna[i]),'tested_windows':path(a.windows[i]),'hyper':path(a.promoters[i*2]),'hypo':path(a.promoters[i*2+1])}
a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(c,indent=2)+'\n')
