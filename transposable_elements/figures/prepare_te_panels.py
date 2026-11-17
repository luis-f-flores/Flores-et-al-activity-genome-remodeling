"""Assemble generated TE measurements for the family and locus plotting programs."""
from pathlib import Path
import argparse,os,shutil

def link(source,destination):
 if not source.is_file():raise FileNotFoundError(source)
 if destination.exists():raise FileExistsError(destination)
 try:os.link(source,destination)
 except OSError:shutil.copyfile(source,destination)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for name in ['context','locus-examples','telescope-run','telescope-de','repeatmasker-de','target-de','family-pca','family-heatmap','copy-proportions','family-null','annotation','output']:
  p.add_argument('--'+name,type=Path,required=True)
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 sources={
  'Fig5_TE_context.tsv.gz':a.context/'te_methylation_context.tsv.gz',
  'Fig5_TE_locus_examples_vector.tsv':a.locus_examples,
  'Fig5_TE_sample_qc.tsv':a.telescope_run/'tables/sample_qc.tsv',
  'Fig5_TE_pca_coordinates.tsv':a.telescope_run/'qc/pca_sample_coordinates.tsv',
  'Fig5_TE_sample_correlation.tsv':a.telescope_run/'qc/sample_correlation_final_count.tsv',
  'Fig5_TE_primary_DE.tsv.gz':a.telescope_de/'results/DE_KCl_vs_Mannitol_final_count.tsv.gz',
  'Fig5_TE_high_confidence_DE.tsv.gz':a.telescope_de/'results/DE_final_conf.tsv.gz',
  'Fig5_TE_subfamily_effects.tsv':a.repeatmasker_de/'tables/repeatmasker_subfamily_DE_KCl_vs_Mannitol.tsv',
  'subfamily_DE_KCl_vs_Mannitol.tsv':a.target_de/'tables/subfamily_DE_KCl_vs_Mannitol.tsv',
  'TE_PCA_all_timepoints_scores.csv':a.family_pca/'TE_PCA_all_timepoints_scores.csv',
  'TE_DE_heatmap_SIGONLY_log2FC_matrix.csv':a.family_heatmap/'TE_DE_heatmap_SIGONLY_log2FC_matrix.csv',
  'TE_DE_heatmap_SIGONLY_stars_matrix.csv':a.family_heatmap/'TE_DE_heatmap_SIGONLY_stars_matrix.csv',
  'Fig5C_copy_proportions.tsv':a.copy_proportions/'Fig5C_copy_proportions.tsv',
  'Fig5C_repeat_family_composition.tsv':a.copy_proportions/'Fig5C_repeat_family_composition.tsv',
  '34_te_family_null_summary.tsv':a.family_null/'tables/34_te_family_null_summary.tsv',
  'SALMONTE_TARGET_ANNOTATION_SUMMARY.tsv':a.annotation/'tables/SALMONTE_TARGET_ANNOTATION_SUMMARY.tsv'}
 for name,source in sources.items():link(source,a.output/name)

if __name__=='__main__':main()
