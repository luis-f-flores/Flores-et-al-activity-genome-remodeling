# DNA methylation

These analyses start from the 21 supplied CpG reports, including baseline, 1 h, 6 h and recovery samples. Reports contain genomic coordinates and methylated/unmethylated read counts.

| Step | Code | Inputs → outputs |
| --- | --- | --- |
| 1. Process CpGs | `cpg_processing/` | CpG reports → coverage summaries, PCA, methylation levels and bigWigs |
| 2. Call changes | `dmr/dmc_1h.R`, `dmr/dmc_6h.R`, `dmr/dmc_recovery.R`, corresponding `dmr_*.R` files and `dmr/merge_dmrs.R` | Per-sample methylation → DMCs, window tests and merged DMRs |
| 3. Annotate intervals | `dmr/build_gtf_features.R`, `dmr/annotate_dmr_features.R`, `dmr/assign_dmr_features.R` | DMRs and gene annotation → genomic-feature overlaps |
| 4. Build sample matrices | `dmr/build_dmr_sample_matrix.py`, `dmr/cpg_coverage_cutoffs.R`, `dmr/count_dmr_regions.cpp`, `dmr/order_dmr_heatmap.py` | Six DMR tables and 21 reports → DMR-by-sample methylation matrix and heatmap order |
| 5. Analyze gene overlap | `gene_overlap/` | DMR/DEG memberships and matching covariates → enrichment tests and GO results |
| 6. Analyze promoters | `promoters/` | Promoter DMRs, RNA tables and eligibility → promoter methylation/RNA cohorts and time-course comparisons |
| 7. Analyze gene bodies | `gene_bodies/build_gene_body_intervals.py`, `gene_bodies/gene_body_methylation.R`, `gene_bodies/gene_body_matching.py`, `gene_bodies/gene_body_metagene_profile.sh` | RefGene GTF/FAI, methylation and RNA → stable gene-body intervals, matched measurements and profiles |
| 8. Relate methylation to RNA | `gene_overlap/prepare_methylation_rna_inputs.py`, `../rna/continuous/` | Coordinate-matched promoter/body measurements and gene counts → TPM associations, OLS regressions and plotting values |
| 9. Draw panels | `figures/`, `loci/` | Generated methylation measurements → Main 4–6 and supplementary panels |

Gene bodies are connected transcript spans from the RefGene annotation. The average profiles use 21,819 loci of at least 2 kb; gene-body methylation tests use the same coordinate construction without a minimum-length filter. RNA matching requires a unique gene assignment with matching name, chromosome, strand and overlapping coordinates.

The root `common/staging.tsv` places the shared helpers where the source-layout scripts expect them. Path and sample helpers are in `methylation/common/`. Use `python run.py plan --dataset wgbs` to inspect the command and input/output contracts. `cpg_processing/align_bismark.sh` records the Read 1 alignment options used before the supplied CpG reports.

The registered preparation stages connect supplied reports and provider RNA tables to the display measurements. `dmr_sample_matrix` generates the DMR union matrix; `wgbs_matrix`, `dmr_heatmap_order`, `dmr_locus_inputs` and `dmr_annotation_summary` feed `wgbs_panel_inputs`. `promoter_dmgs` selects the minimum-q representative before the strict promoter-DMG filter. `gene_body_calls_1`, `gene_body_calls_6` and `gene_body_calls_24` feed the coordinate-matched gene-body/RNA cohorts.

`wgbs_cpg_pca` generates S11B PCA coordinates and its PDF. `wgbs_global_mcg` writes the per-sample and weighted-condition values used in S12A. `gene_body_profiles` takes the generated stable gene-body BED and the three condition-mCG directories to draw S14A. `promoter_locus_inputs` writes the sample methylation, track and RNA values for the Ccdc25 and Serpinb8 examples. The sample bigWigs for these locus panels are the supplied WGBS tracks in `inputs/wgbs/bigwigs/`.
