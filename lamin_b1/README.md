# Lamin B1

This workflow maps the ChIP libraries, measures B1/Input enrichment, fits a shared three-state HMM and relates LAD-associated changes to RNA. Data and generated files stay outside the repository.

| Step | Code | Produces |
| --- | --- | --- |
| Align reads | `chipseq/prepare_chipseq_run.py`, `chipseq/align_chipseq_reads.py` | Duplicate-marked mm10/dm6 BAMs |
| Count fragments | `signal/count_unmasked_bins.py`, `signal/count_sex_chromosome_fragments.py`, `signal/count_b1_libraries.py` | Autosomal and chrX/chrY 5-kb counts with unmasked sequence lengths |
| Normalize B1 | `chipseq/count_spikein_fragments.py`, `signal/prepare_b1_library_totals.py`, `signal/build_b1_signals.py` | Supported 30-kb counts, log2(B1/Input), calibrated signal and linear HMM input |
| Call LADs | `lads/train_shared_hmm.py`, `lads/select_shared_hmm.py`, `lads/summarize_lamin_b1.py` | Shared states, segments, territory, replicate agreement and PCA |
| Map RNA genes | `annotation/build_gene_annotation.py`, `annotation/map_rna_to_refgene.py` | Stable loci and reciprocal-unique RNA/RefGene links |
| Measure gene-span B1 | `transcription/prepare_gene_membership.py`, `transcription/build_lad_gene_measurements.py` | Eligible gene sets, B1 measurements and dependence clusters |
| Test regional B1 | `lads/prepare_shared_lad_comparisons.py`, `lads/regional_limma.R`, `transcription/assemble_pairwise_lad_results.py` | Standardized regional tests and RNA/GO gene links |
| Compare DEGs and controls | `transcription/test_b1_rna_direction.R`, `transcription/prepare_matched_gene_controls.py`, `transcription/matched_b1.R`, `transcription/prepare_lad_overlap_inputs.py`, `transcription/lad_overlap.py` | Directional odds ratios, matched B1 differences and matched DEG/LAD overlaps |
| Summarize nimodipine responses | `transcription/prepare_nimodipine_gene_responses.py`, `transcription/summarize_nimodipine_deg_responses.py` | B1 and RNA responses in the fixed bicuculline DEG sets |
| Profile LADs | `lads/calculate_shared_lad_profiles.py`, `lads/characterize_shared_lads.py` | LAD body/boundary profiles, inside/outside signal and genomic context |
| Compare activity-response genes | `transcription/analyze_activity_gene_lads.py`, `public_lads/compare_activity_gene_lad_maps.py` | Gene-span B1 and LAD-overlap permutations |
| Prepare panel inputs | `signal/prepare_b1_track_panels.py`, `locus_tracks/prepare_locus_tracks.py`, `figures/prepare_b1_figure_data.py` | Genome-wide spread, chr19 and locus arrays, numerical panel tables |
| Draw panels | `figures/` | Main 1–3 and supplementary B1 components |

`chipseq/samples.tsv` identifies treatments, IPs, Inputs and biological replicates. ChIP reads are deposited under GSE354767 and the corresponding RNA data under GSE354768. The gene-expression workflow in [../rna/](../rna/README.md) prepares the RNA tables used here.

Set `data_root`, `work_root` and `output_root` in `config.toml` or a local TOML override. The configured input directory supplies:

```text
lamin_b1/
  fastq/                       # Delivered paired reads listed in raw_reads.tsv
  P2609_gene_fpkm.tsv           # Gene ID, name, biotype and nine per-library FPKM columns
  gene_count.tsv               # Provider counts, including gene_length
references/
  mm10.fa
  mm10.fa.fai
  mm10.refGene.gtf.gz
  genome.gtf.gz                # RNA reference GTF, GRCm38.p6
  org.Mm.eg.sqlite             # org.Mm.eg.db 3.16.0 annotation database
  mm10_blacklist_plus_gaps.bed
  rmsk.txt.gz                 # UCSC mm10 RepeatMasker annotation
```

The annotation builders preserve gene IDs, connected transcript loci, exons and strand-aware TSS coordinates. The RefGene crosswalk requires shared exons, supported identity and reciprocal uniqueness. B1 gene measurements use the full mapped gene span; LAD membership uses its TSS.

From the repository root, inspect the ordered commands and then run the numerical workflow into a new output directory:

```bash
python run.py plan --target b1_pairwise_results --target b1_directional_tests --target matched_b1 --target lad_overlap --target b1_response_summary
python run.py run --target b1_pairwise_results --target b1_directional_tests --target matched_b1 --target lad_overlap --target b1_response_summary --output ../results/b1
```

HMM training uses pomegranate 0.11.1 in the declared legacy environment; set `runtime.legacy_hmm_python` to that interpreter. The launcher runs prerequisites in order. Individual scripts also accept explicit input/output paths, so a stage can be run from existing generated inputs without rebuilding earlier stages. For example:

```bash
python lamin_b1/lads/prepare_shared_lad_comparisons.py --signals /path/b1_signal.npz --decoded /path/selected_decoded.npz --genes /path/gene_coordinates_expression.tsv --rna /path/DEG_classification.tsv --output /path/pairwise
Rscript lamin_b1/lads/regional_limma.R /path/pairwise
python lamin_b1/transcription/assemble_pairwise_lad_results.py --pairwise /path/pairwise --output /path/b1_rna
```

The regional model uses per-IP z scores; tracks, gene-span measurements and signal-spread panels use log2(B1/Input). Shared-HMM calls define the regional analyses, LAD profiles, RNA-heatmap annotations, IAP state distributions and IAP matched controls. The three DEG/LAD overlap comparisons use a common BH correction. Their plotter accepts `--overlap /path/lad_overlap/tables/RESULTS.json`, so annotations come from the calculated results.

After the numerical stages, prepare the display inputs and render the panels:

```bash
python run.py plan --target b1_track_inputs --target b1_locus_inputs --target b1_plot_data
python run.py run --target b1_main1 --target b1_main2 --target b1_main3 --target b1_s04 --target b1_s05 --target b1_s08 --target b1_locus_tracks --target b1_supplement_tracks --output ../results/b1
```

`prepare_continuous_lad_rna.py` writes gene/region links for the continuous RNA workflow. `rna/bic/calculate_baseline_tpm.py` calculates baseline TPM from the provider counts. `prepare_b1_figure_data.py` combines that summary with computed LAD profiles, regional tests, directional tests and GO results. The RNA PCA, expression heatmap and activity scores come from `rna/bic/prepare_bic_expression.py`; the panel renderer accepts multiple `--rna-inputs` directories.

`prepare_locus_tracks.py` calculates the 10-kb B1 tracks and domainograms directly from the fragment-count arrays. RNA curves use condition-mean, 50-bp CPM bigWigs generated by the RNA coverage workflow. `build_smoothed_b1_tracks.py` separately writes the 90-kb smoothed B1 bigWigs from the supported 30-kb counts. Smoothing averages positive linear enrichment over three consecutive supported bins before taking log2; it does not cross unsupported gaps.

Public neural LAD comparisons are in [public_lads/](public_lads/README.md). IAP/B1 analyses are in [../transposable_elements/b1_association/](../transposable_elements/b1_association/).
