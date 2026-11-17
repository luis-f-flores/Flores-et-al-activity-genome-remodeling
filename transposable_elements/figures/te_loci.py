#!/usr/bin/env python3
"""Locus-level TE RNA, methylation and assignment-QC plots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, spearmanr


DATA = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te'
def configure(input_dir: Path | None, output_dir: Path, record_dir: Path) -> None:
    """Set the paths used by the plotting functions."""
    global FIG, RECORD, LOCUS_EXAMPLES, CONTEXT
    global SAMPLE_QC, PCA_COORD, CORR, PRIMARY_DE, CONF_DE, SUBFAMILY_DE
    FIG = output_dir.resolve()
    RECORD = record_dir.resolve()
    root = input_dir.resolve()
    LOCUS_EXAMPLES = root / "Fig5_TE_locus_examples_vector.tsv"
    CONTEXT = root / "Fig5_TE_context.tsv.gz"
    SAMPLE_QC = root / "Fig5_TE_sample_qc.tsv"
    PCA_COORD = root / "Fig5_TE_pca_coordinates.tsv"
    CORR = root / "Fig5_TE_sample_correlation.tsv"
    PRIMARY_DE = root / "Fig5_TE_primary_DE.tsv.gz"
    CONF_DE = root / "Fig5_TE_high_confidence_DE.tsv.gz"
    SUBFAMILY_DE = root / "Fig5_TE_subfamily_effects.tsv"

ORANGE = "#D55E00"
BLUE = "#3B73B9"
GRAY = "#B5B5B5"
DARK = "#333333"

FULL_WIDTH_MM = 170
MAIN_HEIGHT_MM = 235
SUPP_HEIGHT_MM = 178
MM_PER_INCH = 25.4
FULL_WIDTH_IN = FULL_WIDTH_MM / MM_PER_INCH
MAIN_HEIGHT_IN = MAIN_HEIGHT_MM / MM_PER_INCH
SUPP_HEIGHT_IN = SUPP_HEIGHT_MM / MM_PER_INCH

# PCA variance from six Telescope reports using the 6,092-locus filter.
PCA_VARIANCE_PCT = {
    "final_count": (54.957672, 15.200149),
    "final_conf": (53.942038, 15.515733),
}

plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.size": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 9,
})


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def bools(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.lower().eq("true")


def save(fig, stem: str) -> list[Path]:
    pdf = FIG / f"{stem}.pdf"
    png = FIG / f"{stem}.png"
    # Do not use bbox_inches='tight': the page itself is the Illustrator
    # artboard and must remain exactly 170 mm wide.
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    return [pdf, png]


def panel_title(ax, letter: str, title: str, y: float = 1.03,
                title_size: float = 9.0, title_x: float = 0.16) -> None:
    ax.text(0.0, y, letter, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=11, fontweight="bold")
    ax.text(title_x, y, title, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=title_size, fontweight="normal")


def draw_locus_examples(fig, slot, data: pd.DataFrame) -> None:
    'Draw panel G as vector artwork from the compact table.'
    anchor = fig.add_subplot(slot, frame_on=False)
    anchor.set_xticks([]); anchor.set_yticks([]); anchor.patch.set_alpha(0)
    panel_title(anchor, "G", "RNA, methylation and element structure", y=1.025, title_x=0.055)

    # Reserve a row for the panel heading above the locus titles and legends.
    outer = slot.subgridspec(1, 3, width_ratios=[0.065, 0.895, 0.040], wspace=0)
    inner = outer[0, 1].subgridspec(
        4, 5,
        height_ratios=[0.14, 0.71, 0.33, 0.34],
        width_ratios=[1.32, 0.15, 0.70, 0.15, 1.32],
        hspace=0.10,
        wspace=0.18,
    )
    rna_norm = TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)
    mcg_norm = TwoSlopeNorm(vmin=-50, vcenter=0, vmax=50)
    cmap_rna = plt.get_cmap("RdBu_r").copy(); cmap_rna.set_bad("#E5E5E5")
    cmap_mcg = plt.get_cmap("RdBu_r").copy(); cmap_mcg.set_bad("#E5E5E5")

    locus_order = [
        (0, "IAPEZI@chr16:32148444-32149702", "Intronic (Nrros)"),
        (4, "IAPEZI@chr15:39173390-39174662", "Intergenic"),
    ]
    for col, locus_id, title in locus_order:
        z = data[data.selected_locus.eq(locus_id)]
        if z.empty:
            raise RuntimeError(f"Missing compact locus-example records for {locus_id}")
        meta = z.iloc[0]
        start = int(meta.panel_start0); end = int(meta.panel_end0)
        bins = z[z.record_type.eq("mcg_bin")].sort_values("start0")
        components = z[z.record_type.eq("component")].sort_values(["start0", "end0"])
        dmrs = z[z.record_type.eq("dmr")]

        ax = fig.add_subplot(inner[1, col])
        structure = fig.add_subplot(inner[2, col], sharex=ax)
        xedges = np.r_[bins.start0.to_numpy(float), float(bins.end0.iloc[-1])]
        change = np.ma.masked_invalid(bins.mcg_delta_pp.to_numpy(float)[None, :])
        for rr in components.itertuples(index=False):
            x0 = max(float(rr.start0), start); x1 = min(float(rr.end0), end)
            effect = float(rr.rna_log2fc) if np.isfinite(rr.rna_log2fc) else np.nan
            face = "#E5E5E5" if not np.isfinite(effect) else cmap_rna(rna_norm(effect))
            lw = 1.15 if bool(rr.is_selected_component) else 0.55
            edge = "#202020" if bool(rr.is_selected_component) else "white"
            ax.add_patch(Rectangle((x0, 1), x1-x0, 0.58,
                                   facecolor=face, edgecolor=edge, lw=lw))
        ax.pcolormesh(xedges, [0, 1], change, cmap=cmap_mcg, norm=mcg_norm, shading="flat")
        if col == 0:
            ax.set_yticks([1.29, 0.5], ["RNA log2FC", "ΔmCG (%)"], fontsize=9.5)
        else:
            ax.set_yticks([1.29, 0.5], ["", ""])
        ax.set_title(title, loc="left", fontsize=9.5, fontweight="normal", pad=3)
        ax.tick_params(axis="x", labelbottom=False, labelsize=8)
        ax.set_xlim(start, end); ax.set_ylim(0, 1.58)
        ax.spines[["top", "right", "bottom"]].set_visible(False)

        for rr in components.itertuples(index=False):
            label = str(rr.component)
            color = "#9A72AC" if label == "LTR" else "#E9B84F"
            structure.add_patch(Rectangle((rr.start0, 0.13), rr.end0-rr.start0, 0.44,
                                          facecolor=color, edgecolor="white", lw=0.6))
            structure.text((rr.start0+rr.end0)/2, 0.35, label,
                           ha="center", va="center", fontsize=7.2)
        for dr in dmrs.itertuples(index=False):
            structure.add_patch(Rectangle((dr.start0, 0.74), dr.end0-dr.start0, 0.17,
                                          facecolor="#4B97C9", edgecolor="none"))
            structure.text((dr.start0+dr.end0)/2, 0.95, "hypo-DMR",
                           color="#4B97C9", ha="center", va="bottom", fontsize=10.0)
        structure.set_ylim(0, 1.08); structure.set_yticks([]); structure.set_xlim(start, end)
        structure.spines[["top", "right", "left"]].set_visible(False)
        ticks = np.linspace(start, end, 4)
        structure.set_xticks(ticks, [f"{x/1e6:.4f}" for x in ticks], fontsize=7.0)
        structure.set_xlabel(f"{meta.chrom} position (Mb, mm10)", fontsize=9.2, labelpad=3)

    host = data[data.record_type.eq("host_summary")]
    if len(host) != 1:
        raise RuntimeError(f"Expected one host-gene summary, found {len(host)}")
    h = host.iloc[0]
    axg = fig.add_subplot(inner[1:3, 2])
    axg.errorbar([0], [h.host_mean_log2_ratio], yerr=[h.host_sem_log2_ratio],
                 fmt="o", color="#E8742A", ecolor="#E8742A", capsize=5,
                 markersize=7, lw=1.6)
    axg.axhline(0, color="#666666", lw=0.8, linestyle="--")
    axg.set_xlim(-0.48, 0.48); axg.set_xticks([])
    axg.set_ylabel("log2FC", fontsize=9.0, labelpad=1)
    axg.set_title(f"Host-gene mRNA\n<1 FPKM; q={h.host_padj:.3f}",
                  fontsize=9.0, fontweight="normal", pad=4)
    axg.text(0.5, 0.03, f"mean +/- SEM\n{int(h.host_n_pairs)} replicate pairs",
             transform=axg.transAxes, ha="center", va="bottom", fontsize=6.8,
             color="#555555")
    axg.tick_params(axis="y", labelsize=7.5)
    axg.set_ylim(min(-0.35, h.host_mean_log2_ratio-h.host_sem_log2_ratio-0.15),
                 max(1.35, h.host_mean_log2_ratio+h.host_sem_log2_ratio+0.2))
    axg.spines[["top", "right", "bottom"]].set_visible(False)
    axg.grid(axis="y", color="#ECECEC", lw=0.6)

    bars = inner[3, 2].subgridspec(3, 1, height_ratios=[1, 0.95, 1], hspace=0.40)
    cax1 = fig.add_subplot(bars[0, 0]); cax2 = fig.add_subplot(bars[2, 0])
    cb1 = matplotlib.colorbar.ColorbarBase(cax1, cmap=plt.get_cmap("RdBu_r"),
                                           norm=rna_norm, orientation="horizontal")
    cb2 = matplotlib.colorbar.ColorbarBase(cax2, cmap=plt.get_cmap("RdBu_r"),
                                           norm=mcg_norm, orientation="horizontal")
    cax1.set_title("RNA log2FC", fontsize=8.5, pad=2)
    cax2.set_title("ΔmCG (%)", fontsize=8.5, pad=2)
    cax1.tick_params(labelsize=8.0); cax2.tick_params(labelsize=8.0)


def main_components() -> list[Path]:
    pca = read(PCA_COORD)
    primary = read(PRIMARY_DE)[["locus_id", "log2FoldChange_shrunk"]].rename(
        columns={"log2FoldChange_shrunk": "primary"})
    conf = read(CONF_DE)[["locus_id", "log2FoldChange_shrunk"]].rename(
        columns={"log2FoldChange_shrunk": "conf"})
    effects = primary.merge(conf, on="locus_id", how="inner").dropna()
    rho = effects[["primary", "conf"]].corr(method="spearman").iloc[0, 1]
    sub = read(SUBFAMILY_DE).sort_values("log2FoldChange_apeglm").reset_index(drop=True)
    loci = read(CONTEXT)
    loci["iap"] = bools(loci.iap_related)
    top5 = loci.nlargest(5, "rna_log2fc_shrunk").copy()
    locus_examples = read(LOCUS_EXAMPLES)

    fig = plt.figure(figsize=(FULL_WIDTH_IN, MAIN_HEIGHT_IN), facecolor="white")
    gs = fig.add_gridspec(3, 12, height_ratios=[0.92, 0.92, 1.15], hspace=0.34, wspace=0.88,
                          left=0.075, right=0.985, top=0.97, bottom=0.035)

    # A — direct vector workflow.
    ax = fig.add_subplot(gs[0, 0:4]); ax.set_axis_off(); ax.set_box_aspect(1); ax.set_anchor("W")
    pos = ax.get_position()
    ax.set_position([pos.x0 - 0.050, pos.y0, pos.width, pos.height])
    panel_title(ax, "A", "Locus-level TE workflow")
    def workflow_icon(cx, cy, kind):
        if kind == "catalog":
            ax.plot([cx-0.085, cx+0.085], [cy, cy], color="#555555", lw=1.6)
            for dx, color in zip([-0.06, -0.03, 0, 0.03, 0.06],
                                 ["#7B61A8", "#E7B34D", "#6D9DC5", "#999999", ORANGE]):
                ax.add_patch(Rectangle((cx+dx-0.009, cy-0.025), 0.018, 0.05,
                                       facecolor=color, edgecolor="white", lw=0.4))
        elif kind == "samples":
            for i, color in enumerate([ORANGE]*3 + [BLUE]*3):
                xx = cx + (i % 3 - 1) * 0.045
                yy = cy + (0.025 if i < 3 else -0.025)
                ax.scatter(xx, yy, s=38, color=color, edgecolor="white", lw=0.5, zorder=4)
        elif kind == "filter":
            for row in range(2):
                for col in range(6):
                    color = ORANGE if (row == 1 and col in (1, 2, 4)) else "#E4E4E4"
                    ax.add_patch(Rectangle((cx-0.092+col*0.031, cy-0.035+row*0.041), 0.024, 0.031,
                                           facecolor=color, edgecolor="white", lw=0.3))
        elif kind == "pca":
            for dx, dy in [(-.055,.02), (-.025,-.025), (.0,.035)]:
                ax.scatter(cx+dx, cy+dy, s=32, color=BLUE, zorder=4)
            for dx, dy in [(0.035,-.02), (.065,.025), (.09,-.035)]:
                ax.scatter(cx+dx, cy+dy, s=32, color=ORANGE, zorder=4)
        elif kind == "response":
            ax.text(cx-0.025, cy, "↑", ha="center", va="center", color=ORANGE, fontsize=19)
            ax.text(cx+0.035, cy, "↓", ha="center", va="center", color=BLUE, fontsize=19)

    def workflow_box(x, y, w, h, title, detail, kind, title_position="inside", edge="#C9C9C9"):
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.008",
                               facecolor="white", edgecolor=edge, lw=1.0, clip_on=False)
        ax.add_patch(patch)
        workflow_icon(x+w/2, y+h*0.64, kind)
        if title_position == "above":
            title_y, title_va = y + h + 0.010, "bottom"
        elif title_position == "below":
            title_y, title_va = y - 0.008, "top"
        else:
            title_y, title_va = y + h*0.41, "center"
        ax.text(x+w/2, title_y, title, ha="center", va=title_va, fontsize=7.2, linespacing=0.95)
        ax.text(x+w/2, y+h*0.16, detail, ha="center", va="center", fontsize=6.8,
                color="#555555", linespacing=0.95)
    workflow_box(0.00, 0.64, 0.38, 0.21, "RepeatMasker\ncatalog", "395,650 fragments", "catalog", "above")
    workflow_box(0.54, 0.64, 0.38, 0.21, "Any Telescope\nassignment", "191,029 fragments", "samples", "above")
    workflow_box(0.27, 0.35, 0.38, 0.21, "Replicated-support filter", "6,092 loci", "filter", "inside")
    workflow_box(0.00, 0.06, 0.38, 0.20, "PCA", "same 6,092 loci", "pca", "below")
    workflow_box(0.54, 0.06, 0.38, 0.20, "Differential response", "1,049 BH<0.05\n564 higher · 485 lower", "response", "below")
    arrow = dict(arrowstyle="->", color="#888888", lw=1.1, shrinkA=2, shrinkB=2)
    ax.annotate("", (0.54, 0.735), (0.38, 0.735), arrowprops=arrow)
    ax.annotate("", (0.65, 0.455), (0.73, 0.64), arrowprops=arrow)
    ax.annotate("", (0.19, 0.26), (0.38, 0.35), arrowprops=arrow)
    ax.annotate("", (0.73, 0.26), (0.54, 0.35), arrowprops=arrow)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    # B — same-sized direct PCA.
    ax = fig.add_subplot(gs[0, 4:8]); ax.set_box_aspect(1); ax.set_anchor("W")
    pos = ax.get_position()
    ax.set_position([pos.x0 - 0.028, pos.y0, pos.width, pos.height])
    markers = {"A": "o", "B": "s", "C": "^"}
    for row in pca.itertuples():
        color = ORANGE if row.condition == "KCl" else BLUE
        ax.scatter(row.final_count_PC1, row.final_count_PC2, s=62, marker=markers[row.replicate_label],
                   color=color, edgecolor="black", linewidth=0.5, zorder=3)
        ax.annotate(row.replicate_label, (row.final_count_PC1, row.final_count_PC2),
                    xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.axhline(0, color="#DDDDDD", lw=0.8); ax.axvline(0, color="#DDDDDD", lw=0.8)
    pc1, pc2 = PCA_VARIANCE_PCT["final_count"]
    ax.set_xlabel(f"PC1 ({pc1:.1f}%)"); ax.set_ylabel(f"PC2 ({pc2:.1f}%)", labelpad=-6)
    ax.tick_params(axis="y", pad=0)
    panel_title(ax, "B", "TE-locus PCA")
    ax.legend(handles=[Line2D([0],[0], marker="o", color="w", markerfacecolor=BLUE, label="Mannitol", markersize=7),
                       Line2D([0],[0], marker="o", color="w", markerfacecolor=ORANGE, label="KCl", markersize=7)],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), frameon=False,
              ncol=2, fontsize=9, borderaxespad=0.0, columnspacing=1.4)
    style(ax)

    # C — same physical axis size as A and B.
    ax = fig.add_subplot(gs[0, 8:12]); ax.set_box_aspect(1); ax.set_anchor("E")
    ax.scatter(effects.primary, effects.conf, s=5, alpha=0.20, color="#777777", lw=0)
    lim = max(14, np.nanmax(np.abs(effects[["primary", "conf"]].to_numpy())))
    ax.plot([-lim, lim], [-lim, lim], ls="--", lw=1, color="#555555")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Primary log2FC"); ax.set_ylabel("High-confidence log2FC", labelpad=-3)
    panel_title(ax, "C", "Assignment robustness")
    ax.text(0.04, 0.94, f"Spearman ρ={rho:.3f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=8.5)
    style(ax)

    # D — direct ranked RepeatMasker-subfamily estimates.
    ax = fig.add_subplot(gs[1, 0:4])
    colors = {"ERVK": "#C85A5A", "ERV1": "#8C6BB1", "ERVL / MaLR": "#6BAF75",
              "LINE": "#4C78A8", "SINE / Alu": "#E6A23C", "SINE / B2": "#F2C14E"}
    def category(row):
        if row.raw_repeat_class == "LINE": return "LINE"
        if row.raw_repeat_class == "SINE": return f"SINE / {row.raw_repeat_family}"
        if row.raw_repeat_family in ("ERVL", "ERVL-MaLR"): return "ERVL / MaLR"
        if row.raw_repeat_family == "ERV1": return "ERV1"
        if row.raw_repeat_family in ("ERVK", "ERVK?"): return "ERVK"
        raise ValueError(
            f"Unhandled RepeatMasker category: {row.raw_repeat_class}/{row.raw_repeat_family}"
        )
    sub["category"] = [category(r) for r in sub.itertuples()]
    for cat in colors:
        z = sub[sub.category.eq(cat)]
        if len(z): ax.scatter(z.index, z.log2FoldChange_apeglm, s=14, color=colors[cat], label=cat, alpha=0.90)
    ax.axhline(0, color="#777777", lw=0.8)
    labels = ["IAPEz-int", "IAPLTR1_Mm", "IAPEy-int", "IAPA_MM-int", "IAPLTR3-int"]
    label_y = {"IAPEz-int": 2.56, "IAPLTR1_Mm": 2.18, "IAPEy-int": 1.72,
               "IAPA_MM-int": 1.30, "IAPLTR3-int": 0.78}
    for name in labels:
        z = sub[sub.repeatmasker_subfamily.eq(name)]
        if len(z):
            r = z.iloc[0]
            ax.plot([248, r.name - 1.2], [label_y[name], r.log2FoldChange_apeglm],
                    color="#555555", lw=0.8, solid_capstyle="round", zorder=2)
            ax.text(246, label_y[name], name, ha="right", va="center", fontsize=7.2)
    ax.set_xlabel("Subfamily rank"); ax.set_ylabel("Subfamily log2FC", labelpad=-2)
    panel_title(ax, "D", "Subfamily effects")
    ax.legend(frameon=False, fontsize=7.2, ncol=1, loc="upper left",
              handletextpad=0.3, labelspacing=0.28, borderaxespad=0.15)
    style(ax)

    # E — locus-level volcano with the top five labels.
    ax = fig.add_subplot(gs[1, 4:9])
    pos = ax.get_position()
    ax.set_position([pos.x0 + 0.017, pos.y0 + 0.006, pos.width - 0.057, pos.height - 0.012])
    y = -np.log10(np.clip(loci.rna_q_BH.astype(float), 1e-50, 1))
    other = ~loci.iap
    ax.scatter(loci.loc[other, "rna_log2fc_shrunk"], y[other], s=7, color="#C8C8C8", alpha=0.38, lw=0)
    ax.scatter(loci.loc[loci.iap, "rna_log2fc_shrunk"], y[loci.iap], s=12, color=ORANGE, alpha=0.80, lw=0)
    ax.axvline(0, color="#777777", lw=0.8); ax.axhline(-np.log10(0.05), color="#999999", lw=0.8, ls="--")
    key_lines = []
    for number, row in enumerate(top5.itertuples(), start=1):
        q = max(float(row.rna_q_BH), 1e-50)
        yy = -np.log10(q)
        ax.scatter(row.rna_log2fc_shrunk, yy, s=35, color=ORANGE, edgecolor="#7A3300", lw=0.7, zorder=5)
        ax.text(row.rna_log2fc_shrunk, yy, str(number), ha="center", va="center",
                fontsize=6.8, color="white", fontweight="bold", zorder=6)
        key_lines.append(f"{number}  {row.repeat_subfamily}  {row.chrom}:{row.start0/1e6:.2f} Mb")
    ax.set_xlabel("KCl / mannitol shrunken log2FC"); ax.set_ylabel("−log10(q)")
    ax.yaxis.set_label_coords(-0.155, 0.5)
    panel_title(ax, "E", "Locus-level differential expression")
    ax.text(0.02, 0.96, "\n".join(key_lines), transform=ax.transAxes, ha="left", va="top",
            fontsize=6.5, linespacing=1.15,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#DDDDDD", alpha=0.52))
    ax.legend(handles=[Patch(facecolor=ORANGE, label="IAP-related"), Patch(facecolor="#C8C8C8", label="Other TE")],
              frameon=False, fontsize=7.6, loc="upper left", bbox_to_anchor=(0.02, 0.69),
              borderaxespad=0.0, handletextpad=0.45, labelspacing=0.35)
    style(ax)

    draw_locus_examples(fig, gs[2, 0:12], locus_examples)
    return save(fig, "Fig_7_components")


def style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#ECECEC", lw=0.6)


def supplement_components() -> list[Path]:
    """Combine Telescope QC and RNA–WGBS support on one 170-mm artboard."""
    sample = read(SAMPLE_QC)
    pca = read(PCA_COORD)
    corr = read(CORR).rename(columns={"Unnamed: 0": "sample"}).set_index("sample")
    order = ["A6M", "B6M", "C6M", "A6K", "B6K", "C6K"]
    corr = corr.loc[order, order]

    d = read(CONTEXT)
    d["iap"] = bools(d.iap_related)
    d["up"] = bools(d.rna_higher_q_lt_0_05)
    d["hypo"] = bools(d.direct_6hr_hypo_dmr)
    d["hyper"] = bools(d.direct_6hr_hyper_dmr)
    covered = d[(d.shared_covered_positions > 0) & np.isfinite(d.mean_delta_pp_kcl_minus_mannitol)].copy()
    iap_cov = covered[covered.iap].copy()
    focus = covered[covered.iap & covered.up].copy().sort_values("rna_log2fc_shrunk", ascending=False)
    rho_all, p_all = spearmanr(covered.mean_delta_pp_kcl_minus_mannitol, covered.rna_log2fc_shrunk)
    rho_iap, p_iap = spearmanr(iap_cov.mean_delta_pp_kcl_minus_mannitol, iap_cov.rna_log2fc_shrunk)

    fig = plt.figure(figsize=(FULL_WIDTH_IN, SUPP_HEIGHT_IN), facecolor="white")
    outer = fig.add_gridspec(
        2, 1,
        height_ratios=[0.90, 1.0],
        hspace=0.68,
        left=0.105, right=0.95, top=0.95, bottom=0.17,
    )
    top_row = outer[0, 0].subgridspec(1, 3, width_ratios=[0.85, 1.45, 1.10], wspace=0.60)
    bottom_row = outer[1, 0].subgridspec(1, 3, width_ratios=[1.28, 1.30, 0.88], wspace=0.62)

    # A — assignment confidence.
    ax = fig.add_subplot(top_row[0, 0])
    s = sample.set_index("sample").loc[order]
    pct = 100 * s.te_final_conf_sum / s.te_final_count_sum
    colors = [BLUE] * 3 + [ORANGE] * 3
    bars = ax.bar(np.arange(6), pct, color=colors, width=0.72)
    for bar, value in zip(bars, pct):
        ax.text(bar.get_x() + bar.get_width()/2, value - 2.2, f"{value:.1f}",
                ha="center", va="top", rotation=90, fontsize=6.5, color="white")
    ax.set_xticks(np.arange(6), order, rotation=45, ha="right")
    ax.set_ylim(0, 103); ax.set_ylabel("Confident / assigned counts (%)")
    panel_title(ax, "A", "Assignment confidence", y=1.05)
    style(ax)

    # B — high-confidence assignments at the same 6,092 loci.
    ax = fig.add_subplot(top_row[0, 1])
    ax.set_box_aspect(1)
    markers = {"A": "o", "B": "s", "C": "^"}
    for row in pca.itertuples():
        color = ORANGE if row.condition == "KCl" else BLUE
        ax.scatter(row.final_conf_PC1, row.final_conf_PC2, s=36,
                   marker=markers[row.replicate_label], color=color,
                   edgecolor="black", lw=0.45)
        ax.annotate(row.replicate_label, (row.final_conf_PC1, row.final_conf_PC2),
                    xytext=(2, 2), textcoords="offset points", fontsize=6.5)
    ax.axhline(0, color="#DDDDDD", lw=0.7)
    ax.axvline(0, color="#DDDDDD", lw=0.7)
    ax.set_yticks([-25, 0, 25])
    pc1, pc2 = PCA_VARIANCE_PCT["final_conf"]
    ax.set_xlabel(f"PC1 ({pc1:.1f}%)")
    ax.set_ylabel(f"PC2 ({pc2:.1f}%)")
    ax.spines[["top", "right"]].set_visible(False)
    panel_title(ax, "B", "High-confidence PCA", y=1.09)

    # C — sample correlation.
    ax = fig.add_subplot(top_row[0, 2])
    im = ax.imshow(corr.to_numpy(float), cmap="Blues", vmin=0.45, vmax=1.0)
    ax.set_xticks(range(6), order, rotation=45, ha="right")
    ax.set_yticks(range(6), order)
    for i in range(6):
        for j in range(6):
            v = corr.iloc[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=5.4,
                    color="white" if v > 0.72 else DARK)
    panel_title(ax, "C", "Sample correlation\n(Pearson)", y=1.05, title_size=8.2, title_x=0.13)
    cb = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.04)
    cb.set_ticks([0.6, 0.8, 1.0])
    cb.ax.set_title("r", fontsize=9, pad=3)
    cb.ax.tick_params(labelsize=7.5)

    # D — continuous RNA–methylation relationship.
    ax = fig.add_subplot(bottom_row[0, 0])
    ax.scatter(covered.mean_delta_pp_kcl_minus_mannitol, covered.rna_log2fc_shrunk,
               s=5, alpha=0.16, color="#8E8E8E", lw=0,
               label=f"All covered (n={len(covered):,})")
    ax.scatter(iap_cov.mean_delta_pp_kcl_minus_mannitol, iap_cov.rna_log2fc_shrunk,
               s=10, alpha=0.65, color=ORANGE, lw=0,
               label=f"IAP-related (n={len(iap_cov)})")
    ax.axhline(0, color="#777777", lw=0.7); ax.axvline(0, color="#777777", lw=0.7)
    ax.set_xlabel("KCl − Mannitol mean ΔmCG (%)")
    ax.set_ylabel("RNA shrunken log2FC")
    panel_title(ax, "D", "RNA–mCG relationship", y=1.05)
    ax.text(0.03, 0.035, f"All: ρ={rho_all:.3f}, P={p_all:.3f}\nIAP: ρ={rho_iap:.3f}, P={p_iap:.3f}",
            transform=ax.transAxes, fontsize=7.2, va="bottom")
    ax.legend(frameon=False, fontsize=9, loc="upper left", borderaxespad=0.2)
    style(ax)

    # E — complete evidence matrix; every locus has one of two gene-body states.
    ax = fig.add_subplot(bottom_row[0, 1])
    y = np.arange(len(focus))[::-1]
    norm_rna = TwoSlopeNorm(vmin=-10, vcenter=0, vmax=10)
    norm_mcg = TwoSlopeNorm(vmin=-50, vcenter=0, vmax=50)
    for x, val, norm in zip(
        [0, 1, 2],
        [focus.rna_log2fc_shrunk, focus.rna_log2fc_shrunk_final_conf,
         focus.mean_delta_pp_kcl_minus_mannitol],
        [norm_rna, norm_rna, norm_mcg],
    ):
        ax.scatter(np.repeat(x, len(focus)), y, c=val, cmap="RdBu_r", norm=norm,
                   marker="s", s=18, lw=0)
    dmr_color = np.where(focus.hypo, "#2C7FB8", np.where(focus.hyper, "#D95F0E", "#D9D9D9"))
    inside = focus.gene_context.eq("genic")
    context_color = np.where(inside, "#7B3294", "#E6E6E6")
    ax.scatter(np.repeat(3, len(focus)), y, c=dmr_color, marker="s", s=18, lw=0)
    ax.scatter(np.repeat(4, len(focus)), y, c=context_color, marker="s", s=18,
               edgecolor="#BDBDBD", linewidth=0.15)
    ax.set_xticks(range(5), ["Primary", "High-conf.", "ΔmCG", "DMR", "Gene body"],
                  rotation=45, ha="right", fontsize=8.2)
    ax.set_yticks([]); ax.set_xlim(-0.6, 4.6); ax.set_ylim(-1, len(focus) + 7)
    panel_title(ax, "E", "61-locus evidence matrix", y=1.05)
    ax.spines[:].set_visible(False)

    # Compact continuous scales for the two RNA columns and the methylation
    # column.  Keep the DMR/gene-body key categorical.
    rna_cax = ax.inset_axes([0.02, 0.925, 0.40, 0.028])
    mcg_cax = ax.inset_axes([0.55, 0.925, 0.40, 0.028])
    rna_cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm_rna, cmap="RdBu_r"),
                          cax=rna_cax, orientation="horizontal")
    mcg_cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm_mcg, cmap="RdBu_r"),
                          cax=mcg_cax, orientation="horizontal")
    rna_cb.set_ticks([-10, 0, 10]); mcg_cb.set_ticks([-50, 0, 50])
    rna_cb.ax.tick_params(labelsize=7, length=2, pad=1)
    mcg_cb.ax.tick_params(labelsize=7, length=2, pad=1)
    rna_cax.set_title("RNA log2FC", fontsize=7.4, pad=2)
    mcg_cax.set_title("ΔmCG (%)", fontsize=7.4, pad=2)
    matrix_legend = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#2C7FB8", label="hypo-DMR", markersize=6),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#D95F0E", label="hyper-DMR", markersize=6),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#7B3294", label="inside gene body", markersize=6),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#E6E6E6", markeredgecolor="#BDBDBD",
               label="outside gene body", markersize=6),
    ]
    ax.legend(handles=matrix_legend, frameon=False, fontsize=9, loc="upper center",
              bbox_to_anchor=(0.5, -0.22), ncol=2, handletextpad=0.3, columnspacing=0.7)

    # F — concise two-state gene-body comparison.
    ax = fig.add_subplot(bottom_row[0, 2])
    counts = [(int((inside & focus.hypo).sum()), int(inside.sum())),
              (int((~inside & focus.hypo).sum()), int((~inside).sum()))]
    vals = [100 * a / b for a, b in counts]
    bars = ax.bar([0, 1], vals, color=["#7B3294", "#BDBDBD"], width=0.62)
    for bar, (a, b), v in zip(bars, counts, vals):
        ax.text(bar.get_x()+bar.get_width()/2, v+0.8, f"{v:.1f}%\n{a}/{b}",
                ha="center", fontsize=8.2)
    _, fisher_p = fisher_exact([[counts[0][0], counts[0][1]-counts[0][0]],
                                [counts[1][0], counts[1][1]-counts[1][0]]])
    ax.set_xticks([0, 1], ["Inside", "Outside\ngene bodies"])
    ax.set_ylabel("RNA-up IAPs with hypo-DMR (%)")
    ax.set_ylim(0, max(vals) + 8)
    panel_title(ax, "F", "Gene-body\ncontext", y=1.05)
    ax.text(0.98, 0.94, f"Fisher P={fisher_p:.3f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=8.2)
    style(ax)

    return save(fig, "Fig_S17_components")


def record(path: Path) -> dict:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": h.hexdigest()}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DATA,
                        help="Directory containing the TE numerical tables.")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Directory for the two PDF/PNG figure pairs.")
    parser.add_argument("--record-dir", type=Path,
                        help="Optional secondary output directory.")
    return parser.parse_args()


def main():
    args = parse_args()
    configure(args.input_dir, args.output_dir, args.record_dir or args.output_dir)
    FIG.mkdir(parents=True, exist_ok=True)
    RECORD.mkdir(parents=True, exist_ok=True)
    inputs = [LOCUS_EXAMPLES, CONTEXT, SAMPLE_QC, PCA_COORD, CORR,
              PRIMARY_DE, CONF_DE, SUBFAMILY_DE]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    outputs = []
    outputs += main_components()
    outputs += supplement_components()

    print('Generated TE locus and assignment-QC plots.')


if __name__ == "__main__":
    main()
