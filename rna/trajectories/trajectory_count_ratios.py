#!/usr/bin/env python3



from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)


import pandas as pd
import numpy as np
import os


input_dir = os.path.expanduser(gb_path('rna/K_clustering_Neuron_RNA_seq_11-13-24/'))
output_dir = os.path.join(input_dir, "filtered_results/")
os.makedirs(output_dir, exist_ok=True)

log2fc_threshold = 0.37
padj_threshold = 0.01


timepoints = {
    "Hour_1_KClvsHour_1_Mann_deg.xlsx": ("1hr", "Hour_1_KCl", "Hour_1_Mann"),
    "Hour_6_KClvsHour_6_Mann_deg.xlsx": ("6hr", "Hour_6_KCl", "Hour_6_Mann"),
    "Hour_24_KClvsHour_24_Mann_deg_all_numeric.xlsx": ("24hr", "Hour_24_KCl", "Hour_24_Mann")
}


all_degs = set()
gene_info = {}
raw_dfs = {}


for file_name, (label, kcl_col, mann_col) in timepoints.items():
    file_path = os.path.join(input_dir, file_name)
    df = pd.read_excel(file_path)
    raw_dfs[label] = df

    
    passing = df[(df["padj"] <= padj_threshold) & (df["log2FoldChange"].abs() >= log2fc_threshold)]
    all_degs.update(passing["gene_id"].tolist())

    for _, row in passing.iterrows():
        gene_info[row["gene_id"]] = row["gene_name"]


final_frames = []

for label, df in raw_dfs.items():
    _, kcl_col, mann_col = [v for k, v in timepoints.items() if v[0] == label][0]

    
    subset = df[df["gene_id"].isin(all_degs)].copy()

    
    subset[kcl_col] = subset[kcl_col].fillna(0)
    subset[mann_col] = subset[mann_col].fillna(0)
    subset[f"log2FC_{label}"] = np.log2((subset[kcl_col] + 1) / (subset[mann_col] + 1))

    trimmed = subset[["gene_id", "gene_name", f"log2FC_{label}"]].drop_duplicates().set_index("gene_id")
    final_frames.append(trimmed)


combined = pd.concat(final_frames, axis=1, join="outer")


combined["gene_name"] = combined.filter(like="gene_name").bfill(axis=1).iloc[:, 0]
combined = combined.loc[:, ~combined.columns.str.contains("gene_name") | (combined.columns == "gene_name")]


log2fc_cols = [col for col in combined.columns if col.startswith("log2FC_")]
combined = combined[~(combined[log2fc_cols].fillna(0).eq(0).all(axis=1))]


output_file = os.path.join(output_dir, "combined_log2FC_FPKM_GlobalDEG.csv")
combined.to_csv(output_file)
print(f"Normalized-count trajectory matrix saved: {output_file}")
