# Telescope locus expression

The paired-read workflow uses Cutadapt, Bowtie2 `--very-sensitive-local -k 500 --score-min L,0,1.6`, query-name sorting and Telescope 1.0.3. The final EM fit uses epsilon `1e-6`, up to 100,000 iterations, theta prior 200,000, pi prior 0, reassign mode `exclude` and confidence 0.9. Libraries are fitted separately.

Prepare an external reference directory containing `mm10.fa.fai`, the RepeatMasker-fragment `te_fragments.gtf` and the six Bowtie2 index files under `bowtie2/mm10`. The TE annotation programs are in `../annotation/`. Use a separate environment containing Telescope 1.0.3 for the annotation cache and quantification:

```bash
/path/to/telescope/bin/python telescope_annotation_cache.py --build-cache /path/to/references/telescope_annotation.pickle --gtf /path/to/references/te_fragments.gtf
python telescope_full_library.py --stage prepare --samples /path/to/samples.tsv --reference-dir /path/to/references --telescope-python /path/to/telescope/bin/python --run-dir /path/to/telescope_run
python telescope_full_library.py --stage run --run-dir /path/to/telescope_run
```

The sample table has `sample`, `fastq1`, `fastq2`, `condition` and `replicate_label` columns. FASTQ paths may be relative to the table. The six-hour comparison uses conditions `Mannitol` and `KCl`, with replicate labels `A`, `B` and `C`. The runner trims paired reads, aligns chunks, merges all candidate alignments and fits one complete library at a time. Disk-backed loading preserves the Telescope mapping order. Completed stages can be resumed from the same run directory.

The output contains each library's Telescope report, compressed reports under `source_reports/`, sample metadata and `tables/sample_qc.tsv`. For the six-hour KCl comparison:

```bash
python prepare_telescope_de_inputs.py --source-package /path/to/telescope_run --output /path/to/telescope_de
Rscript telescope_differential_expression.R /path/to/telescope_de
python summarize_telescope_de.py --results-root /path/to/telescope_de
python telescope_sample_quality.py --package-root /path/to/telescope_run
```

The DE program uses the mapped-fragment size factors, final counts and high-confidence counts in the generated matrices. `telescope_sample_quality.py` calculates the sample PCA and correlations from these reports. RepeatMasker subfamily and TE-target summaries use the same generated DE input directory:

```bash
Rscript telescope_subfamily_expression.R /path/to/telescope_de /path/to/TE_annotation_master_mm10.tsv.gz /path/to/repeatmasker_expression
Rscript telescope_target_expression.R /path/to/telescope_de /path/to/TE_meta.tsv /path/to/SalmonTE_out_6h/A6M/quant.sf /path/to/SalmonTE_DE_6h/TE_DE_results.csv /path/to/target_expression
```
