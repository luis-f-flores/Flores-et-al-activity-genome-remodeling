"""Compare activity-response gene coverage across the six public LAD maps."""
from pathlib import Path
import argparse
import importlib.util
import numpy as np
import pandas as pd

MAPS = [('E13 Dorsal NPC', 'Dorsal_NPCs_merged_lads.bed'),
        ('E13 Ventral NPC', 'Ventral_NPCs_merged_lads.bed'),
        ('30wk Hippocampus', 'LaminB1_Hippocampus_LAD_peaks_30wk_mm10.bed'),
        ('E13 Liver', 'Liver_merged_lads.bed'), ('3T3 MEF', 'MEF_merged_lads.bed'),
        ('E13 Lung LaminB', 'combined_E13_Lung_LaminB.bed')]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['coordinates', 'gene-lists', 'maps', 'output']:
        p.add_argument('--' + key, type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).resolve().parents[1] / 'transcription/analyze_activity_gene_lads.py'
    spec = importlib.util.spec_from_file_location('activity_lads', source)
    helper = importlib.util.module_from_spec(spec); spec.loader.exec_module(helper)
    genes = pd.read_csv(a.coordinates, sep='\t', names=['chr', 's', 'e', 'gene'])
    if not genes.gene.is_unique or not (genes.e > genes.s).all():
        raise ValueError('Gene coordinates must be unique, positive-length intervals.')
    membership, indices = helper.memberships(genes, a.gene_lists)
    membership.to_csv(a.output / 'activity_gene_membership.tsv', sep='\t', index=False)
    fractions = []
    for label, filename in MAPS:
        intervals = pd.read_csv(a.maps / filename, sep='\t', header=None).iloc[:, :3]
        intervals.columns = ['chr', 's', 'e']
        lookup = helper.interval_lookup(intervals)
        fractions.append([helper.overlap(r.chr, r.s, r.e, lookup) / (r.e-r.s) for r in genes.itertuples()])
    flags = np.array(fractions).T >= .75
    result, nulls = helper.permutations(flags, indices, [x[0] for x in MAPS])
    result.to_csv(a.output / 'public_LAD_permutation.tsv', sep='\t', index=False)
    np.savez_compressed(a.output / 'public_LAD_permutation_nulls.npz', **nulls)
    for j, (label, _) in enumerate(MAPS): genes[label + '_LAD_fraction'] = fractions[j]
    genes.to_csv(a.output / 'public_LAD_gene_coverage.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
