"""Resolve CpG report filenames and merged DMR paths for methylation workflows."""
from pathlib import Path
import argparse,csv
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--reports',type=Path,required=True);p.add_argument('--dmr-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
def write(name,columns,rows):
 with (a.output/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader();w.writerows(rows)
rows=[]
for sid in [r+'0' for r in 'ABC']+[f'{r}{t}{c}' for t in [1,6,24] for c in ['M','K'] for r in 'ABC']:
 paths=list(a.reports.glob(sid+'_*.CpG_report.txt'))+list(a.reports.glob(sid+'.CpG_report.txt'))
 if len(paths)!=1:raise ValueError(f'Expected one CpG report for {sid}; found {len(paths)}')
 rows.append({'sample':sid,'report':str(paths[0].resolve())})
write('samples.tsv',['sample','report'],rows);write('controls.tsv',['sample','report'],[r for r in rows if r['sample'].endswith('M')])
dmrs=[]
for t in [1,6,24]:
 for dr in ['hyper','hypo']:
  path=a.dmr_root/f'DMR_{t}M_vs_{t}K_multiwindow/DMR_win_500bp/Merged_q05_meth5/{t}M_vs_{t}K_Merged_DMRs_cov4_{dr}.txt'
  dmrs.append(dict(timepoint=t,direction=dr,path=str(path.resolve())))
write('merged_dmrs.tsv',['timepoint','direction','path'],dmrs)
