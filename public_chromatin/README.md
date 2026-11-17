# Public chromatin and baseline methylation

These programs compare public postnatal cortical DNMT3A/H3K9me2 measurements with TTX LADs and baseline mCG.

| Step | Code | Inputs → outputs |
| --- | --- | --- |
| 1. Count reads | `processing/count_public_libraries.py`, `processing/count_public_chromatin_reads.py` | Deposited read coordinates and chromosome sizes → ChIP/Input counts on the original grid |
| 2. Project to LADs | `processing/project_chromatin_to_lads.py` | Counts, mm9-to-mm10 chain, genomic masks and TTX LADs → LAD enrichment and body profiles |
| 3. Aggregate CpGs | `processing/prepare_b1_context.py`, `processing/extract_public_mcg.py`, `processing/aggregate_baseline_mcg.py`, `processing/aggregate_mcg_dnmt3a_windows.py` | Generated B1 states, public and culture CpG reports → TTX LAD tables and methylation counts in B1 and DNMT3A windows |
| 4. Calculate comparisons | `analysis/chromatin_mcg_windows.py`, `analysis/chromatin_context.py` | Projected chromatin and CpG counts → rank groups, body profiles and correlations |
| 5. Draw panels | `figures/chromatin_context.py` | Calculated summaries → Supplementary Figure 10 panels |

The coordinate projection uses the UCSC chain's aligned blocks. Count normalization, unmapped positions, eligibility and averaging are implemented in the processing scripts. DNMT3A uses all three deposited ChIP/Input pairs. H3K9me2 uses the deposited ChIP library and both H3 controls. Baseline culture mCG uses the three baseline WGBS samples.

The read-coordinate filenames and sample roles are listed in `processing/samples.tsv`. Place those files from GSE104578 in `inputs/public_chromatin/reads/`. The public CpG files are `GSM2757470_Cortex_BS_1wk.txt.gz` and `GSM2757471_Cortex_BS_2wk.txt.gz`, in `inputs/public_chromatin/methylation/`. The reference directory contains `mm9.chrom.sizes`, `mm9ToMm10.over.chain.gz`, `mm10.fa` and its FASTA index. Culture CpG reports use the WGBS working layout described in `methylation/README.md`.

```bash
python run.py plan --figure S10
python run.py run --target chromatin_context_plot --data-root /path/to/inputs --work-root /path/to/work --output /path/to/results
```

The workflow builds the B1 grid and TTX LADs from the ChIP measurements and selected HMM, counts the public reads, aggregates CpGs and generates the plotted summaries.
