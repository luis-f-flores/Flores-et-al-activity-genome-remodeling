"""Calculate the genic DMRâ€“DEG overlap null using the supplied covariates and strata."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'genic_overlap');p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    results=[]
    for tp in [1,6,24]:
        membership=pd.read_csv(a.inputs/f'{tp}h_membership.tsv',sep='\t')
        strata=pd.read_csv(a.inputs/f'{tp}h_strata.tsv',sep='\t')
        np.testing.assert_array_equal(np.sort(membership.stratum.unique()),strata.stratum)
        rng=np.random.default_rng(20261020+tp);draws=np.zeros(10000,dtype=int)
        for r in strata.itertuples():
            m=membership[membership.stratum==r.stratum];assert len(m)==r.N and m.DEG.sum()==r.n and m.DMR.sum()==r.DMR_genes
            if r.n:draws+=rng.hypergeometric(r.DMR_genes,r.N-r.DMR_genes,r.n,size=10000)
        observed=int((membership.DEG&membership.DMR).sum());expected=float(membership.expected_inclusion@membership.DMR.to_numpy(int));P=float((1+(draws>=observed).sum())/10001)
        results.append(dict(timepoint=tp,observed=observed,expected=expected,P=P))
        pd.DataFrame({'overlap':draws}).to_csv(a.output/f'{tp}h_null.tsv',sep='\t',index=False)
    (a.output/'RESULTS.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results))

if __name__=='__main__':main()
