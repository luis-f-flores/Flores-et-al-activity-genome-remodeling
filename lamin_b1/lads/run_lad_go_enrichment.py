'Run regional GO enrichment with explicit file paths.'
from pathlib import Path
import argparse,os,subprocess
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--membership',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--rscript',default=os.environ.get('RSCRIPT','Rscript'))
a=p.parse_args();source=Path(__file__).with_name('lad_go_enrichment.R')
a.output.mkdir(parents=True,exist_ok=True)
subprocess.run([a.rscript,'--vanilla',str(source),str(a.membership),str(a.output),str(source)],check=True)
