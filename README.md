# Neuronal activity reshapes genome–lamina interactions and the DNA methylome alongside distinct transcriptional responses

Analysis code for Lamin B1 ChIP-seq, RNA-seq, WGBS, transposable elements and DamID in primary cortical neurons. Each assay directory follows the measurements from input processing to statistical analysis and figure panels.

## Workflows

| Directory | Analysis |
| --- | --- |
| [lamin_b1/](lamin_b1/README.md) | ChIP alignment and enrichment, LAD calling, transcriptional relationships and locus tracks |
| [methylation/](methylation/README.md) | CpG reports, DMCs/DMRs, promoter and gene-body methylation |
| [rna/](rna/README.md) | RNA abundance, temporal clusters and KCl–bicuculline comparisons |
| [transposable_elements/](transposable_elements/README.md) | Family and locus expression, IAP methylation and B1 association |
| [damid/](damid/README.md) | GATC counting, continuous enrichment, read-depth matching and ranked windows |
| [public_chromatin/](public_chromatin/README.md) | Postnatal DNMT3A/H3K9me2 and baseline methylation comparisons |

## Run the code

Use Python 3.12. Install Python packages from `requirements.txt`; R packages and other tools are listed in [environment/](environment/). Sequencing-processing commands run on Linux, including WSL.

```bash
python -m pip install -r requirements.txt
python run.py list --dataset chipseq
python run.py plan --figure Fig_1
python run.py check --target kcl_tpm --data-root /path/to/inputs
python run.py run --target kcl_tpm --data-root /path/to/inputs --output /path/to/results
```

`config.toml` sets input, working and output directories and software commands. `--config` accepts a local TOML override. `--data-root`, `--work-root` and `--output` override individual paths. The default directories are siblings of this repository.

`plan` shows a target's exact inputs, outputs and command. `check` verifies its prerequisites. `run` executes the selected calculations in dependency order. Processing stages marked `source` are run with the commands in their assay guides; their generated inputs must exist before selecting downstream plots. For the original sequencing-processing layout, `python run.py prepare --work-root /path/to/work` places the source scripts and configuration files in the required working directories. The assay guides describe the processing order and input files.

The [figure workflow table](figure_workflows.tsv) connects each manuscript figure to its analysis and plotting code. [workflows.tsv](workflows.tsv) defines runnable targets and their dependencies. Figure programs write numerical tables and computational panels; [common/staging.tsv](common/staging.tsv) defines the processing-code layout.

## Data

Keep input data and generated results in the configured external directories. Sequencing data are deposited under:

| Assay | GEO Series |
| --- | --- |
| DamID | [GSE329354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329354) |
| KCl RNA-seq | [GSE329355](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329355) |
| WGBS | [GSE329357](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329357) |
| Lamin B1 ChIP-seq | [GSE354767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354767) |
| Bicuculline RNA-seq | [GSE354768](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354768) |

WGBS starts from supplied Bismark CpG reports. Gene RNA analyses start from the provider's count, gene-length and differential-expression tables. ChIP and TE processing start from sequencing reads; the DamID counter accepts aligned reads and a GATC annotation. Reference assemblies are mouse mm10/GRCm38 and Drosophila dm6. Sample identifiers containing `24` denote recovery.

## Licence and citation

Original code is available under the [MIT licence](LICENSE). Third-party DamID code retains its licence and attribution in [THIRD_PARTY.md](THIRD_PARTY.md). Citation metadata are in [CITATION.cff](CITATION.cff).
