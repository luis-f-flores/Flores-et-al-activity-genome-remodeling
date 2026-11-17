# Public neural LAD maps

The workflow uses the supplied UCSC mm10 RefSeq gene-coordinate table, the rPRG/dPRG/SRG gene lists and six public LAD maps. Keep these reference files under the configured external `data_root/public_lads` directory.

| Step | Program | Output |
| --- | --- | --- |
| Prepare gene intervals | `extract_gtf_coordinates.py` | `cleaned_gene_coordinates.bed` |
| Summarize map territory | `summarize_lad_maps.py` | LAD sequence coverage |
| Compare activity-response genes | `compare_activity_gene_lad_maps.py` | Gene overlap fractions and permutation results |
| Draw comparisons | `plot_public_lad_enrichment.py` | Public-LAD enrichment panels |

The coordinate input is the five-column table `mm10_refseq_USCS_gene_coordinates_only.gtf`: chromosome, strand, gene start, gene end and symbol. The builder preserves its coordinate convention and row order while removing noncanonical chromosome names. The supplied gene lists are `rPRG_gene_list.txt`, `dPRG_gene_list.txt` and `SRG_gene_list.txt`; Fam46a and Gltscr2 resolve to Tent5a and Nop53.

The six BED inputs are dorsal and ventral E13 NPCs, adult hippocampus, E13 liver, MEFs and E13 lung. Their filenames are declared in `compare_activity_gene_lad_maps.py` and `workflows.tsv`. A gene overlaps a LAD when at least 75% of its complete annotated span is covered by the union of LAD intervals. The comparison uses 10,000 sampled gene sets from the same gene-coordinate universe.

```bash
python run.py plan --target public_lad_plot
python run.py run --target public_lad_plot --output ../results/public_lads
```

Individual stages also accept paths directly:

```bash
python lamin_b1/public_lads/extract_gtf_coordinates.py --coordinates /path/mm10_refseq_USCS_gene_coordinates_only.gtf --output /path/genes.bed
python lamin_b1/public_lads/compare_activity_gene_lad_maps.py --coordinates /path/genes.bed --gene-lists /path/gene_lists --maps /path/lad_beds --output /path/enrichment
python lamin_b1/public_lads/plot_public_lad_enrichment.py --input /path/enrichment/public_LAD_permutation.tsv --fonts /path/fonts --output /path/panels
```

The same interval and permutation functions are used by `../transcription/analyze_activity_gene_lads.py` for the experimental B1 maps.
