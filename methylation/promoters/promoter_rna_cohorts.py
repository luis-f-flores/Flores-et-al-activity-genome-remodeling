"""Match selected promoter DMRs to provider RNA records using stable gene identifiers."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale

def model_rows(df):
    sub=df.dropna(subset=['meth_diff','log2FC_fpkm']).copy()
    assert np.isfinite(sub[['meth_diff','log2FC_fpkm']].to_numpy()).all()
    sub['meth_diff_scaled']=minmax_scale(sub['meth_diff'].to_numpy(),feature_range=(-.5,.5))
    return sub

def match_promoter_genes(old, rna, overrides, tp):
    output, changes = [], []
    numeric = ['mean_KCl', 'mean_Mann', 'log2FC_fpkm', 'FC_fpkm', 'padj']
    for symbol, previous in old.groupby('gene_name', sort=False):
        entry = previous.iloc[0].to_dict()
        override = overrides[overrides['dmg_symbol'].eq(symbol)]
        if len(override):
            assert len(override) == 1
            resolved_id = override.iloc[0]['rna_gene_id']
            candidates = rna[rna['gene_id'].eq(resolved_id)]
            reason = override.iloc[0]['reason']
        else:
            candidates = rna[rna['gene_name'].eq(symbol)]
            if len(candidates) > 1:
                raise ValueError(f'Unresolved multiple RNA IDs for {symbol} at {tp}hr')
            reason = 'exact unique symbol'
        if len(candidates) > 1:
            raise ValueError(f'Nonunique stable RNA ID for {symbol}')
        for col in ['gene_id'] + numeric:
            entry[col] = np.nan
        entry['rna_gene_name'] = None
        entry['join_method'] = reason
        if len(candidates):
            row = candidates.iloc[0]
            k = np.mean([row[f'{rep}{tp}K_fpkm'] for rep in 'ABC'])
            m = np.mean([row[f'{rep}{tp}M_fpkm'] for rep in 'ABC'])
            fc = (k+1)/(m+1)
            fresh = dict(gene_id=row['gene_id'], padj=row['padj'], mean_KCl=k,
                         mean_Mann=m, FC_fpkm=fc, log2FC_fpkm=np.log2(fc))
            entry.update(fresh)
            entry['rna_gene_name'] = row['gene_name']
        output.append(entry)
    result = pd.DataFrame(output)
    assert result['gene_name'].is_unique
    assert result['gene_id'].dropna().is_unique
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--promoter-files',type=Path,nargs=6,required=True,help='Hyper/hypo promoter-DMG tables for 1 h, 6 h and recovery.')
    ap.add_argument('--rna',type=Path,nargs=3,required=True)
    ap.add_argument('--overrides',type=Path,default=Path(__file__).with_name('promoter_gene_id_overrides.tsv'))
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    overrides=pd.read_csv(a.overrides,sep='\t');counts=[]
    for i,tp in enumerate([1,6,24]):
        d=pd.concat([pd.read_csv(p,sep='\t') for p in a.promoter_files[2*i:2*i+2]],ignore_index=True)
        assert d.gene_name.is_unique
        assert (d.qvalue<=.01).all() and (d.meth_diff.abs()>=25).all()
        rna=pd.read_csv(a.rna[i],sep='\t')
        joined=match_promoter_genes(d,rna,overrides,tp);model=model_rows(joined)
        joined.to_csv(a.output/f'20_DMG_gene_matched_{tp}hr.tsv',sep='\t',index=False)
        model.to_csv(a.output/f'20_DMG_model_input_{tp}hr_scaled.tsv',sep='\t',index=False)
        sig=rna[(rna.padj<=.01)&(rna.log2FoldChange.abs()>=np.log2(1.3))]
        counts.append(dict(timepoint=f'{tp}hr',all_DMGs=joined.gene_name.nunique(),DEGs=sig.gene_name.nunique(),corrected_overlap=len(set(joined.gene_id.dropna())&set(sig.gene_id))))
    pd.DataFrame(counts).to_csv(a.output/'promoter_RNA_counts.tsv',sep='\t',index=False)

if __name__=='__main__':main()
