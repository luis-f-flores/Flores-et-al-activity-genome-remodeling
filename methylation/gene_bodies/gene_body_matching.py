#!/usr/bin/env python3
"""Match methylated gene bodies to RNA genes by name, chromosome, strand and overlap."""
import argparse
import gzip
import re
from pathlib import Path
import numpy as np
import pandas as pd


def chromosome(value):
    value = str(value).removeprefix('chr')
    return 'GL456233.1' if value == 'X_GL456233_random' else value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rna-gtf', type=Path, required=True)
    parser.add_argument('--gene-bed', type=Path, required=True)
    parser.add_argument('--gene-body', type=Path, nargs=3, required=True)
    parser.add_argument('--rna', type=Path, nargs=3, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    reference = {}
    opener = gzip.open if args.rna_gtf.suffix == '.gz' else open
    with opener(args.rna_gtf, 'rt') as stream:
        for line in stream:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] != 'gene':
                continue
            attrs = dict(re.findall(r'(\w+) "([^\"]*)"', fields[8]))
            gene_id = attrs['gene_id']
            assert gene_id not in reference
            reference[gene_id] = (attrs.get('gene_name'), chromosome(fields[0]),
                                  int(fields[3]), int(fields[4]), fields[6])
    bed = pd.read_csv(args.gene_bed, sep='\t', header=None,
                     names=['chr','start0','end','locus_id','score','strand','gene_symbol'])
    bed_keys = set(zip(bed.chr, bed.start0+1, bed.end, bed.strand, bed.locus_id, bed.gene_symbol))
    candidates, rejected, accounting = [], [], []
    for tp, body_path, rna_path in zip(('1hr','6hr','24hr'), args.gene_body, args.rna):
        body = pd.read_csv(body_path, float_precision='round_trip')
        rna = pd.read_csv(rna_path, sep='\t', float_precision='round_trip')
        assert rna.gene_id.is_unique
        selected = body[(body['meth.diff'].abs()>=5) & (body.qvalue<=.05)].copy()
        assert selected.locus_id.is_unique
        assert all(k in bed_keys for k in zip(selected.chr,selected.start,selected.end,
                   selected.strand,selected.locus_id,selected.gene_symbol))
        columns = ['gene_id','gene_name','gene_chr','gene_start','gene_end','gene_strand','log2FoldChange','padj']
        joined = selected.merge(rna[columns], left_on='gene_symbol', right_on='gene_name')
        for row in joined.drop_duplicates('gene_id').itertuples():
            assert (row.gene_name,chromosome(row.gene_chr),row.gene_start,row.gene_end,row.gene_strand) == reference[row.gene_id]
        valid = ((joined.chr.map(chromosome)==joined.gene_chr.map(chromosome)) &
                 (joined.strand==joined.gene_strand) & (joined.start<=joined.gene_end) &
                 (joined.end>=joined.gene_start))
        rejected.append(joined[~valid].assign(reason='incompatible_RNA_coordinates'))
        candidates.append(joined[valid])
        accounting.append(selected[['timepoint','locus_id','gene_symbol','chr','start','end','strand']])
    candidates = pd.concat(candidates, ignore_index=True)
    mapping = candidates[['locus_id','gene_id']].drop_duplicates()
    ambiguous_loci = set(mapping.loc[mapping.locus_id.duplicated(False),'locus_id'])
    ambiguous_rna = set(mapping.loc[mapping.gene_id.duplicated(False),'gene_id'])
    ambiguous = candidates.locus_id.isin(ambiguous_loci) | candidates.gene_id.isin(ambiguous_rna)
    rejected.append(candidates[ambiguous].assign(reason='nonunique_locus_RNA_assignment'))
    matched = candidates[~ambiguous].rename(columns={'gene_id':'gene_id_rna'}).copy()
    matched['body_locus_id'] = matched.locus_id
    crosswalk = matched[['body_locus_id','gene_id_rna']].drop_duplicates()
    assert crosswalk.body_locus_id.is_unique and crosswalk.gene_id_rna.is_unique
    ledger = pd.concat(accounting, ignore_index=True)
    retained = set(zip(matched.timepoint,matched.locus_id))
    ledger['RNA_match_status'] = ['retained' if (r.timepoint,r.locus_id) in retained
                                  else 'no_unique_coordinate_match' for r in ledger.itertuples()]
    args.output.mkdir(parents=True, exist_ok=True)
    matched.to_csv(args.output/'64_displayed_gene_locus_RNA_matches.csv', index=False)
    crosswalk.to_csv(args.output/'64_unique_locus_to_RNA_ID_crosswalk.csv', index=False)
    pd.concat(rejected,ignore_index=True).to_csv(args.output/'64_unmatched_RNA_candidates.csv', index=False)
    ledger.to_csv(args.output/'64_all_DMG_match_accounting.csv', index=False)
    print('Matched gene bodies:',matched.groupby('timepoint').size().to_dict())


if __name__ == '__main__':
    main()
