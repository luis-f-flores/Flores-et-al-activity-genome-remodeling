"""Write file paths for TE-family shuffles and IAP genomic matching."""
from pathlib import Path
import argparse
import json
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['context','rna-de','dmr-root','dmr-manifest','gene-layers','mcg','bismap','reports','chrom-sizes','exclusions','te-meta','family-overlap','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    def absolute(path):return str(Path(path).resolve())
    table=pd.read_csv(args.dmr_manifest,sep='\t')
    dmrs={}
    for t in [1,6,24]:
        rows=table[table.timepoint.eq(t)]
        if len(rows)!=2 or set(rows.direction)!={'hyper','hypo'}:raise ValueError(f'Invalid DMR manifest at {t} h')
        dmrs[f'{t}hr']={r.direction:absolute(r.path) for r in rows.itertuples()}
        dmrs[f'{t}hr']['assessed']=absolute(args.dmr_root/f'DMR_{t}M_vs_{t}K_multiwindow/DMR_win_500bp/DMRs_all_unfiltered_win_500bp.txt')
    family=dict(chrom_sizes=absolute(args.chrom_sizes),exclusions=absolute(args.exclusions),te_meta=absolute(args.te_meta),family_overlap=absolute(args.family_overlap),dmrs=dmrs)
    (args.output/'family_null.json').write_text(json.dumps(family,indent=2)+'\n')
    iap=json.loads(Path(__file__).with_name('iap_settings.json').read_text())
    iap.update(context=absolute(args.context),rna=absolute(args.rna_de),dmrs=dmrs['6hr']['hypo'],tested_windows=dmrs['6hr']['assessed'],
               genes=absolute(args.gene_layers/'genes.bed'),exons=absolute(args.gene_layers/'exons.bed'),promoters=absolute(args.gene_layers/'promoters_-1kb_+100bp.bed'),
               kcl_bw=absolute(args.mcg/'16_genomewide_KCl_6hr_mean_pct_meth.bw'),mann_bw=absolute(args.mcg/'16_genomewide_Mannitol_6hr_mean_pct_meth.bw'),
               bismap=absolute(args.bismap),raw_reports=absolute(args.reports))
    (args.output/'iap_matching.json').write_text(json.dumps(iap,indent=2)+'\n')


if __name__=='__main__':main()
