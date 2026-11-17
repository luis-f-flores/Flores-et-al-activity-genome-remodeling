# Transposable elements

Family and locus measurements follow separate quantification routes, then connect to methylation and Lamin B1. Inputs and results are stored outside this repository.

| Stage | Entry points | Input and output |
| --- | --- | --- |
| Reference catalog | `annotation/build_te_catalog_mm10.sh`, `catalog_copy_identifiers.sh`, `annotate_te_catalog_mm10.py` | RepeatMasker and SalmonTE reference sequences produce fragment coordinates, identifiers and the annotation master |
| Family expression | `family_expression/salmonte_quantify.sh`, `family_de_*.R`, `family_pca.R`, `te_differential_heatmap.R` | RNA reads produce SalmonTE counts, differential expression, PCA and heatmap values |
| Bic family comparison | `family_expression/bic_salmonte_quantification.py`, `salmonte_cross_paradigm.R` | Salmon quantifications and read totals produce IAP replicate CPM, family changes and DESeq2 contrasts |
| Locus expression | `telescope/` | Paired RNA reads produce Telescope fragment counts, sample PCA and locus/subfamily DESeq2 results; see its [guide](telescope/README.md) |
| Methylation context | `methylation/te_methylation_context.py` | Locus expression, mCG bigWigs, DMRs and gene layers produce TE methylation/RNA/context tables |
| Locus examples | `methylation/prepare_locus_examples.py` | TE context, DMR overlaps, annotation, mCG bigWigs and gene RNA tables produce the two IAP example panels |
| Family overlap | `methylation/te_dmr_overlap.R`, `family_dmr_null.py` | The TE catalog and merged DMRs produce family overlap counts and a chromosome/length-matched null in the assessed methylation space |
| IAP enrichment | `methylation/run_iap_matching.py` | IAP fragments, CpG support and genomic references produce matched positions and fragment/grouped hypo-DMR enrichment |
| IAP/B1 association | `b1_association/` | TE coordinates and generated B1 signals produce LAD overlap and matched B1 measurements |
| Plotting | `figures/prepare_te_panels.py`, `te_loci.py`, `te_families.py` | Generated family, locus and methylation tables produce Main 7 and supplementary components |

The IAP workflow constructs candidate positions, measures CpG support, completes the matching pools and calculates enrichment. `te_matching_configs` resolves the generated context, DMR, gene-layer and condition-mCG inputs together with the six CpG reports and genome references. The numerical matching settings are in `methylation/iap_settings.json`.

```bash
python transposable_elements/methylation/run_iap_matching.py --settings /path/to/iap_settings.json --output /path/to/results/iap
Rscript transposable_elements/family_expression/salmonte_cross_paradigm.R /path/to/KCl_RNA /path/to/Bic_quant /path/to/results/family_comparison
```

The family comparison reads each library's `quant.sf` and `aux_info/meta_info.json`. The KCl directory contains `SalmonTE_out_1h`, `SalmonTE_out_6h`, `SalmonTE_out_24h` and the corresponding `SalmonTE_DE_*` directories. The Bic directory contains the nine sample quantification directories.

The IAP matching tables are written under `results/`, and its candidate pools under `candidate_pools/`. `family_dmr_null.py` uses 1,000 seeded permutations; its configuration names the chromosome sizes, exclusion BED, TE metadata, observed family table and tested/merged DMR tables. Use `python run.py plan --dataset te` for the registered commands and dependencies.

`te_catalog_fragments` builds the coordinate catalog and `te_catalog` annotates its fragments. `telescope_de_inputs` consumes the reports produced by the raw-read workflow, then `telescope_de`, `telescope_repeatmasker_de`, `telescope_target_de` and `telescope_sample_profiles` generate the expression and sample tables. `te_context`, `te_locus_examples`, `te_family_dmr_null`, `te_copy_proportions` and `te_family_heatmap` provide the remaining inputs to `te_panel_inputs`.
