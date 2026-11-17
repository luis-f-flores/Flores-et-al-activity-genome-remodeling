# Gene expression

The RNA analyses use provider count matrices, gene lengths and DESeq2 tables for the KCl time course and the TTX/bicuculline/nimodipine experiment.

| Step | Code | Inputs → outputs |
| --- | --- | --- |
| 1. Build temporal trajectories | `trajectories/trajectory_count_ratios.py` | Normalized count tables → temporal expression ratios |
| 2. Cluster and annotate | `trajectories/cluster_trajectories.R`, `trajectories/cluster_go.R` | Trajectories → cluster assignments and GO enrichment |
| 3. Calculate TPM displays | `expression/order_trajectory_heatmap.R`, `expression/kcl_tpm.py` | Counts, gene lengths, clusters and activity-gene membership → original-order TPM row Z scores |
| 4. Compare stimulation paradigms | `expression/prepare_rna_comparison.py`, `expression/rna_comparison.py` | Provider DE tables and counts → jointly testable genes, expression changes and Spearman correlations |
| 5. Relate RNA to chromatin | `continuous/continuous_rna.py`, `continuous/prepare_continuous_panels.py` | Counts, methylation cohorts and gene-linked LADs → correlations, OLS slopes, R² and regression curves |
| 6. Draw panels | `figures/` | Calculated TPM matrices and cross-paradigm comparisons → figure panels |

Continuous expression change is log2[(mean target TPM + 1)/(mean reference TPM + 1)]. DESeq2 calls retain the provider's count-based analysis. The B1/RNA spatial calculation uses 10,000 within-chromosome permutations of 5-Mb blocks.

Use `python run.py plan --dataset rna` for the target commands. [trajectories/PROVIDER.md](trajectories/PROVIDER.md) describes the provider-table starting point.

## Bicuculline and nimodipine

`bic/prepare_bic_expression.py` calculates size factors, transcriptome PCA, activity-class scores and DEG summaries from the provider tables. Its heatmap stage clusters the bicuculline DEGs and adds TTX LAD annotations from the shared-HMM calls. Place `gene_count.xls`, `BicvsTTX_deg.xls`, `Bic_NimovsBic_deg.xls` and `Bic_NimovsTTX_deg.xls` in `work/rna_package/inputs/P2609-01_gene_tables/`; the rPRG, dPRG and SRG lists go in `work/rna_package/genesets/`.

```bash
python run.py plan --target bic_expression --target bic_heatmap
python run.py run --target bic_expression --output ../results
```

`bic/build_rna_coverage.py` creates locus-display tracks from the nine supplied, coordinate-sorted BAMs and their indexes in `inputs/rna/bic/bams/`. Sample names are `Basal_A`–`Basal_C`, `Bic_A`–`Bic_C` and `Bic_nimo_A`–`Bic_nimo_C`. Coverage uses 50-bp bins and primary mapped read ends per million retained read ends, combining strands and retaining duplicate-marked reads. Each condition track averages its three separately normalized libraries. The B1 locus workflow consumes these generated tracks.

```bash
python run.py plan --target bic_rna_coverage
python run.py run --target bic_rna_coverage --output ../results
```

`expression_inputs`, `activity_membership`, `kcl_gene_order` and `kcl_tpm_inputs` prepare the count, gene-length and ordering tables consumed by `kcl_tpm`. `rna_testability` uses the provider KCl and Bic-versus-TTX tables; `rna_comparison` calculates the expression changes from generated counts. `continuous_inputs` combines generated promoter/body cohorts with the gene-linked B1 measurements before `continuous_rna` calculates its associations.
