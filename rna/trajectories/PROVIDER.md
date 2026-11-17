# Gene RNA processing

Novogene supplied the KCl and Bic gene count/expression and differential-expression tables used in this study. The recorded workflow used FastQC, Trim Galore, HISAT2 alignment to mm10 and DESeq2. The analyses in this repository begin with those supplied tables. Novogene did not provide the executable pipeline used to generate them.

The supplied count matrices and gene lengths support the TPM calculations. DEG selection uses the supplied DESeq2 results and the thresholds specified in Methods. The original clustering uses normalized-count ratios, row standardization and Ward.D2; the TPM display retains those gene memberships and ordering. Large provider reports and raw libraries remain in the sequence-data deposit.
