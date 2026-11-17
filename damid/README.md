# DamID

The DamID workflow compares Dam–Lamin B1 with Dam-only controls at 6 h and recovery, with two biological replicates per treatment and construct.

| Step | Code | Inputs → outputs |
| --- | --- | --- |
| 1. Count assigned reads | `processing/count_gatc_batch.sh`, `processing/gatc_counting_mm10.R` | R1 alignments, GATC annotation and chromosome sizes → per-feature and 400-kb counts |
| 2. Aggregate support | `processing/prepare_window_counts.py` | GATC and 400-kb counts → count arrays, positive-site support and library totals |
| 3. Calculate continuous measurements | `continuous/measure_continuous_signal.py` | Count arrays and library totals → 2-Mb enrichment, read support, depth-matched SD and ranked changes |
| 4. Select plotted windows | `continuous/select_rank_examples.py` | Ranked windows and gene annotation → annotated rank-band callouts and locus summaries |
| 5. Draw panels | `figures/damid_continuous.py` | Continuous measurements → read-support, SD, rank and locus panels |

Place the 16 R1 BAMs listed in `processing/samples.example.tsv` in `inputs/damid/bams/`, named by sample ID. Place the mm10 GATC GFF and chromosome sizes in `inputs/references/`. The original R counter assigns reads to GATC features and the 400-kb grid. `legacy_gene_layers` builds the RefGene intervals used to annotate the windows from the reference GTF.

The 2-Mb sliding windows advance by 400 kb. The depth-matched SD uses the specified random seed and 1,000 subsamples; rank plots use full-depth continuous signals. Window tables and plotted values are generated from the assigned counts.

```bash
python run.py plan --target damid_continuous_plot
python run.py run --target damid_continuous_plot --data-root /path/to/inputs --output /path/to/results
```

The third-party GATC counter retains its licence in `processing/LICENSE`.
