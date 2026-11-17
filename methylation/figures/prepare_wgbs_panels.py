"""Assemble generated DMR matrices, annotations and locus values for WGBS panels."""
from pathlib import Path
import argparse,os,shutil

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['matrix','pca','ordering','loci','annotation','gene-layers','output']:p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 sources={
  'DMR_sample_matrix.npz':a.matrix/'DMR_sample_matrix.npz',
  'PCA_scores_PC1_PC2_PC3.tsv':a.pca/'DMR_PCA_scores.tsv',
  'PCA_variance.tsv':a.pca/'DMR_PCA_variance.tsv',
  'heatmap_cluster_order.tsv.gz':a.ordering/'heatmap_cluster_order.tsv.gz',
  'three_example_plot_values.tsv':a.loci/'three_example_plot_values.tsv',
  'three_example_DMRs.tsv':a.loci/'three_example_DMRs.tsv',
  'original_DMR_annotation_pies.json':a.annotation/'original_DMR_annotation_pies.json',
  'exons.bed':a.gene_layers/'exons.bed'}
 for name,source in sources.items():
  if not source.is_file():raise FileNotFoundError(source)
  destination=a.output/name
  if destination.exists():raise FileExistsError(destination)
  try:os.link(source,destination)
  except OSError:shutil.copyfile(source,destination)

if __name__=='__main__':main()
