"""Match the published activity-response gene lists to provider RNA identifiers."""
from pathlib import Path
import argparse
import pandas as pd

ALIASES={'5430416O09Rik':'Spaar','9430020K01Rik':'Jcad','Fam46a':'Tent5a','Gltscr2':'Nop53'}

def prepare(rna,lists,output):
    genes=pd.read_csv(rna,sep='\t',usecols=['gene_id','gene_name'])
    rows=[]
    for group,path in zip(['rapidPRGs','delayedPRGs','SRGs'],lists):
        for symbol in dict.fromkeys(s.strip() for s in Path(path).read_text().splitlines() if s.strip()):
            original={v:k for k,v in ALIASES.items()}.get(symbol,symbol)
            name=ALIASES.get(symbol,symbol);q=genes[genes.gene_name.eq(name)]
            if len(q)>1:raise ValueError(f'Multiple RNA IDs for activity gene {symbol}')
            if len(q):rows.append(dict(original_symbol=original,gene_name=name,gene_id=q.gene_id.iloc[0],gene_set=group,alias_restored=original!=name))
    d=pd.DataFrame(rows).sort_values('original_symbol',kind='stable').reset_index(drop=True)
    if not d.gene_id.is_unique:raise ValueError('Activity classes contain repeated RNA identifiers.')
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(output,sep='\t',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--rna',type=Path,required=True)
    p.add_argument('--lists',type=Path,nargs=3,required=True,help='rPRG, dPRG and SRG reference gene lists.')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();prepare(a.rna,a.lists,a.output)
