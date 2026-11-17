#!/usr/bin/env python3
'Draw the 170-mm SalmonTE family figure from eight numerical tables.\n\nThe caller supplies a flat input directory and an output directory.'

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


INPUT_FILENAMES = {
    "pca": "TE_PCA_all_timepoints_scores.csv",
    "heatmap_values": "TE_DE_heatmap_SIGONLY_log2FC_matrix.csv",
    "heatmap_stars": "TE_DE_heatmap_SIGONLY_stars_matrix.csv",
    "copy_proportions": "Fig5C_copy_proportions.tsv",
    "copy_composition": "Fig5C_repeat_family_composition.tsv",
    "dmr_null": "34_te_family_null_summary.tsv",
    "telescope_de": "subfamily_DE_KCl_vs_Mannitol.tsv",
    "annotation": "SALMONTE_TARGET_ANNOTATION_SUMMARY.tsv",
}

OUTDIR = Path(".")
PCA_SCORES = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["pca"])
HEATMAP_VALUES = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["heatmap_values"])
HEATMAP_STARS = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["heatmap_stars"])
COPY_PROPORTIONS = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["copy_proportions"])
COPY_COMPOSITION = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["copy_composition"])
R34_TABLE = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["dmr_null"])
R42_TABLE = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["telescope_de"])
ANNOTATION_TABLE = Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'te' / Path(INPUT_FILENAMES["annotation"])

PAGE_W_MM = 170.0
PAGE_H_MM = 220.0
PT = 1.0
BASE_FONT = 9.0
TITLE_FONT = 9.5
PANEL_FONT = 12.0

NAVY = "#172f3b"
GRID = "#e6eaec"
CONTROL = "#f27f79"
KCL = "#2bb8c1"
RED = "#cf4f4f"
ORANGE = "#d96a00"
BLUE = "#4f6fa8"
GRAY = "#87959d"
DARK_GRAY = "#52666f"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def axes_mm(fig: plt.Figure, x: float, y: float, width: float, height: float):
    return fig.add_axes(
        [x / PAGE_W_MM, y / PAGE_H_MM, width / PAGE_W_MM, height / PAGE_H_MM]
    )


def panel_letter(fig: plt.Figure, letter: str, x: float, y: float) -> None:
    fig.text(
        x / PAGE_W_MM,
        y / PAGE_H_MM,
        letter,
        ha="left",
        va="top",
        fontsize=PANEL_FONT,
        fontweight="bold",
        color=NAVY,
    )


def style_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#44545d")
    ax.spines["bottom"].set_color("#44545d")
    ax.tick_params(labelsize=BASE_FONT, length=3, color="#44545d")


def plot_panel_a(fig: plt.Figure) -> None:
    ax = axes_mm(fig, 10, 160, 52, 50)
    data = pd.read_csv(PCA_SCORES)
    colors = {"control": CONTROL, "KCl": KCL}
    markers = {1: "o", 6: "^", 24: "s"}
    for (condition, time), group in data.groupby(["condition", "time"], sort=False):
        ax.scatter(
            group["PC1"],
            group["PC2"],
            s=34,
            marker=markers[int(time)],
            color=colors[str(condition)],
            edgecolor="white",
            linewidth=0.5,
            zorder=3,
        )
    ax.axhline(0, color="#aab4b9", lw=0.8)
    ax.axvline(0, color="#aab4b9", lw=0.8)
    ax.grid(color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    ax.set_xlim(-2.25, 4.75)
    ax.set_ylim(-3.9, 1.9)
    ax.set_xlabel("PC1 (24.9%)", fontsize=BASE_FONT)
    ax.set_ylabel("PC2 (14.7%)", fontsize=BASE_FONT)
    ax.set_title("TE PCA (1 h + 6 h + 24 h)", fontsize=TITLE_FONT, fontweight="normal", pad=4)
    style_axes(ax)

    condition_handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=5.5,
               markerfacecolor=colors[name], markeredgecolor="white", label=name)
        for name in ("control", "KCl")
    ]
    time_handles = [
        Line2D([], [], marker=markers[t], linestyle="none", markersize=5.5,
               markerfacecolor="#222222", markeredgecolor="#222222", label=str(t))
        for t in (1, 6, 24)
    ]
    legend1 = ax.legend(
        handles=condition_handles,
        title="Condition",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.02),
        frameon=False,
        borderaxespad=0,
        handletextpad=0.35,
        labelspacing=0.25,
        fontsize=BASE_FONT,
        title_fontsize=BASE_FONT,
    )
    ax.add_artist(legend1)
    ax.legend(
        handles=time_handles,
        title="Time (h)",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.57),
        frameon=False,
        borderaxespad=0,
        handletextpad=0.35,
        labelspacing=0.25,
        fontsize=BASE_FONT,
        title_fontsize=BASE_FONT,
    )
    panel_letter(fig, "A", 3, 216)


def plot_panel_b(fig: plt.Figure) -> None:
    values = pd.read_csv(HEATMAP_VALUES)
    stars = pd.read_csv(HEATMAP_STARS, keep_default_na=False)
    if not values[["family", "type"]].equals(stars[["family", "type"]]):
        raise RuntimeError("Heatmap value/star row identities differ")
    effect_columns = ["1h KCl vs Mann", "6h KCl vs Mann", "24h KCl vs Mann"]
    # Place the last CSV row at the top.
    values = values.iloc[::-1].reset_index(drop=True)
    stars = stars.iloc[::-1].reset_index(drop=True)
    matrix = values[effect_columns].to_numpy(float)

    ax = axes_mm(fig, 116, 103, 48, 107)
    cmap = LinearSegmentedColormap.from_list(
        "te_diverging", ["#5a70b8", "#f7f7f7", "#e34a33"], N=256
    )
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-3, vmax=3, interpolation="nearest")
    ax.set_yticks(np.arange(len(values)), values["family"], fontsize=8.0)
    ax.set_xticks(np.arange(3), ["1 h", "6 h", "24 h"], fontsize=BASE_FONT)
    ax.set_xlabel("KCl vs mannitol", fontsize=BASE_FONT, labelpad=3)
    ax.set_title("TE expression by subfamily", fontsize=TITLE_FONT, fontweight="normal", pad=4)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row in range(len(values)):
        for col, column in enumerate(effect_columns):
            label = str(stars.loc[row, column])
            if label:
                ax.text(col, row, label, ha="center", va="center", fontsize=8.0, color="#111111")

    cax = axes_mm(fig, 129, 92, 24, 3.5)
    cb = fig.colorbar(image, cax=cax, orientation="horizontal", ticks=[-3, -2, -1, 0, 1, 2, 3])
    cb.ax.tick_params(labelsize=BASE_FONT, length=2, pad=1)
    cb.outline.set_visible(False)
    cb.set_label("log2FC", fontsize=BASE_FONT, labelpad=2)
    panel_letter(fig, "B", 74, 216)


def plot_panel_c(fig: plt.Figure) -> dict[str, object]:
    data = pd.read_csv(COPY_PROPORTIONS, sep="\t")
    composition = pd.read_csv(COPY_COMPOSITION, sep="\t")
    expected = {
        "Higher in KCl": 23957,
        "Lower in KCl": 4740,
        "Not significant": 1135111,
    }
    observed = dict(zip(data["status"], data["copies"].astype(int)))
    if observed != expected:
        raise RuntimeError(f"Panel C verified totals changed: {observed}")

    expected_higher_composition = {
        "ERVK": 20047,
        "ERVL": 2685,
        "L1": 840,
        "ERV1": 385,
    }
    higher = composition.loc[composition["status"].eq("Higher in KCl")]
    observed_higher_composition = dict(
        zip(higher["repeat_family"], higher["copies"].astype(int))
    )
    if observed_higher_composition != expected_higher_composition:
        raise RuntimeError(
            "Panel C verified higher-in-KCl clade composition changed: "
            f"{observed_higher_composition}"
        )

    ax = axes_mm(fig, 7, 93, 70, 49)
    ax.set_xlim(-1.25, 1.75)
    ax.set_ylim(-1.05, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    counts = np.array(
        [expected["Not significant"], expected["Lower in KCl"], expected["Higher in KCl"]],
        dtype=float,
    )
    colors = ["#d9d9d9", BLUE, "#ef1b23"]
    wedges, _ = ax.pie(
        counts,
        colors=colors,
        startangle=90,
        counterclock=True,
        radius=0.58,
        center=(-0.60, 0.03),
        wedgeprops={"linewidth": 0.6, "edgecolor": "#555555"},
    )
    ax.text(-0.60, -0.34, "97.53%", ha="center", va="center", fontsize=BASE_FONT)
    ax.annotate(
        "2.06%",
        xy=(-0.53, 0.61),
        xytext=(-0.34, 0.84),
        ha="center",
        va="center",
        fontsize=BASE_FONT,
        arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#444444"},
    )
    ax.annotate(
        "0.41%",
        xy=(-0.61, 0.61),
        xytext=(-0.92, 0.77),
        ha="center",
        va="center",
        fontsize=BASE_FONT,
        arrowprops={"arrowstyle": "-", "lw": 0.6, "color": "#444444"},
    )

    # Show the RepeatMasker family/clade composition of the higher-in-KCl group.
    clade_order = ["ERVK", "ERVL", "L1", "ERV1"]
    clade_colors = {
        "ERVK": "#ef1b23",
        "ERVL": "#f27f83",
        "L1": "#f6cdcf",
        "ERV1": "#f6d4ad",
    }
    bar_x = 0.19
    bar_w = 0.15
    bar_bottom = -0.36
    bar_height = 0.98
    cumulative = bar_bottom
    for clade in reversed(clade_order):
        fraction = expected_higher_composition[clade] / expected["Higher in KCl"]
        segment_height = bar_height * fraction
        ax.bar(
            bar_x,
            segment_height,
            width=bar_w,
            bottom=cumulative,
            color=clade_colors[clade],
            edgecolor="#555555",
            linewidth=0.45,
        )
        cumulative += segment_height
    ax.plot([-0.53, bar_x - bar_w / 2], [0.61, bar_bottom + bar_height], color="#444444", lw=0.6)
    ax.plot([-0.61, bar_x - bar_w / 2], [0.61, bar_bottom], color="#444444", lw=0.6)

    ax.text(
        0.48,
        0.72,
        "Higher-in-KCl\nclades",
        ha="left",
        va="center",
        fontsize=BASE_FONT,
    )

    for index, clade in enumerate(clade_order):
        y_position = 0.41 - index * 0.20
        ax.add_patch(
            plt.Rectangle(
                (0.48, y_position - 0.045),
                0.09,
                0.09,
                facecolor=clade_colors[clade],
                edgecolor="#555555",
                linewidth=0.45,
            )
        )
        ax.text(0.61, y_position, clade, ha="left", va="center", fontsize=BASE_FONT)

    direction_items = [
        ("#d9d9d9", "Not changed"),
        (BLUE, "Lower"),
        ("#ef1b23", "Higher"),
    ]
    direction_x = [-1.13, 0.05, 0.72]
    for x_position, (color, label) in zip(direction_x, direction_items):
        ax.scatter(
            [x_position + 0.05],
            [-0.71],
            marker="s",
            s=36,
            facecolor=color,
            edgecolor="#555555",
            linewidth=0.45,
            zorder=4,
        )
        ax.text(x_position + 0.13, -0.71, label, ha="left", va="center", fontsize=BASE_FONT)
    ax.text(-0.48, -0.93, "total = 1,163,808 copies", ha="center", va="center", fontsize=BASE_FONT)
    ax.set_title(
        "6 h clade-level TE\ndifferential expression",
        fontsize=TITLE_FONT,
        fontweight="normal",
        pad=3,
    )
    panel_letter(fig, "C", 3, 149)
    return {
        "catalog_copies": int(sum(expected.values())),
        "higher_in_kcl": int(expected["Higher in KCl"]),
        "lower_in_kcl": int(expected["Lower in KCl"]),
        "not_significant": int(expected["Not significant"]),
        "higher_in_kcl_clade_copies": expected_higher_composition,
    }


def plot_panel_d(fig: plt.Figure) -> dict[str, object]:
    summary = pd.read_csv(R34_TABLE, sep="\t")
    data = summary.loc[summary["null_model"].eq("callable")].copy()
    data = data.sort_values("observed_pct", ascending=True).reset_index(drop=True)
    data["above_95"] = data["observed_pct"] > data["null_q975_pct"]
    data["below_95"] = data["observed_pct"] < data["null_q025_pct"]
    if len(data) != 20 or int(data["above_95"].sum()) != 10 or int(data["below_95"].sum()) != 0:
        raise RuntimeError("Unexpected callable-background classification")

    ax = axes_mm(fig, 30, 11, 55, 63)
    y = np.arange(len(data))
    ax.hlines(y, data["null_q025_pct"], data["null_q975_pct"], color="#9aa6ac", lw=1.4, zorder=1)
    ax.scatter(data["null_mean_pct"], y, s=18, color=DARK_GRAY, zorder=2)
    ax.scatter(data["observed_pct"], y, s=20, color=RED, zorder=3)
    for row_index, row in data.loc[data["above_95"]].iterrows():
        ax.text(
            row["observed_pct"] + 0.38,
            row_index,
            "*",
            ha="left",
            va="center",
            fontsize=BASE_FONT,
            fontweight="bold",
            color=NAVY,
            zorder=4,
        )
    ax.set_yticks(y, data["TE_family"], fontsize=BASE_FONT)
    ax.set_xlabel("Family copies overlapping >=1 DMR (%)", fontsize=BASE_FONT)
    ax.set_title(
        "TE-DMR overlap relative to\nWGBS-callable background",
        fontsize=TITLE_FONT,
        fontweight="normal",
        pad=4,
    )
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    # Reserve a compact key area above the top family so the legend does not
    # obscure any observed/null interval.
    ax.set_ylim(-0.8, 23.0)
    style_axes(ax)
    handles = [
        Line2D([], [], marker="o", linestyle="-", color="#9aa6ac", markerfacecolor=DARK_GRAY,
               markeredgecolor=DARK_GRAY, markersize=4, label="Null (95%)"),
        Line2D([], [], marker="o", linestyle="none", color=RED, markerfacecolor=RED,
               markersize=4, label="Observed"),
        Line2D([], [], marker="$*$", linestyle="none", color=NAVY, markersize=6,
               label="Above 95%"),
    ]
    ax.legend(
        handles=handles,
        frameon=False,
        loc="upper left",
        fontsize=BASE_FONT,
        handlelength=1.2,
        handletextpad=0.45,
        labelspacing=0.25,
        borderaxespad=0.15,
    )
    panel_letter(fig, "D", 3, 86)
    return {
        "families": int(len(data)),
        "above_central_95": int(data["above_95"].sum()),
        "below_central_95": int(data["below_95"].sum()),
    }


def read_iap_lookup() -> dict[str, bool]:
    lookup: dict[str, bool] = {}
    with ANNOTATION_TABLE.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lookup[row["salmonte_target"].lower()] = int(row["n_iap_loci"]) > 0
    return lookup


def plot_panel_e(fig: plt.Figure) -> dict[str, float | int]:
    data = pd.read_csv(R42_TABLE, sep="\t")
    data = data.loc[
        np.isfinite(data["salmonte_6h_log2FoldChange"])
        & np.isfinite(data["log2FoldChange_apeglm"])
    ].copy()
    if len(data) != 187:
        raise RuntimeError(f"Expected 187 matched effects; found {len(data)}")
    lookup = read_iap_lookup()
    data["iap_related"] = data["repeat_subfamily"].str.lower().map(lookup).fillna(False)
    x = data["salmonte_6h_log2FoldChange"].to_numpy(float)
    y = data["log2FoldChange_apeglm"].to_numpy(float)
    iap = data["iap_related"].to_numpy(bool)
    rho, p_value = spearmanr(x, y)
    iap_rho, iap_p_value = spearmanr(x[iap], y[iap])
    if not math.isclose(float(rho), 0.4825028399, abs_tol=5e-10):
        raise RuntimeError("All-target Spearman value changed")
    if not math.isclose(float(iap_rho), 0.8877192982, abs_tol=5e-10):
        raise RuntimeError("IAP-target Spearman value changed")

    ax = axes_mm(fig, 108, 14, 57, 59.4)
    limit = 3.2
    ax.axhline(0, color="#9ca9b0", lw=0.7, linestyle=(0, (3, 3)), zorder=0)
    ax.axvline(0, color="#9ca9b0", lw=0.7, linestyle=(0, (3, 3)), zorder=0)
    ax.plot([-3, 3], [-3, 3], color="#c5cdd1", lw=0.8, linestyle=(0, (5, 4)), zorder=0)
    ax.scatter(x[~iap], y[~iap], s=18, color=GRAY, alpha=0.65, linewidth=0,
               label=f"Other TE target (n={int((~iap).sum())})", zorder=2)
    ax.scatter(x[iap], y[iap], s=25, color=ORANGE, alpha=0.95, linewidth=0.35,
               edgecolor="white", label=f"IAP-related target (n={int(iap.sum())})", zorder=3)

    target_rows = data.set_index("repeat_subfamily")
    label_specs = {
        "IAPEZI": (5, 6, "left"),
        "IAPLTR1_MM": (5, -18, "left"),
        "RLTR13B3": (-5, 6, "right"),
    }
    for target, (dx, dy, align) in label_specs.items():
        row = target_rows.loc[target]
        color = "#b85600" if bool(row["iap_related"]) else "#465a65"
        ax.annotate(
            target,
            (float(row["salmonte_6h_log2FoldChange"]), float(row["log2FoldChange_apeglm"])),
            xytext=(dx, dy),
            textcoords="offset points",
            ha=align,
            va="center",
            fontsize=BASE_FONT,
            color=color,
            arrowprops={"arrowstyle": "-", "color": color, "lw": 0.55},
            zorder=4,
        )
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("SalmonTE log2 fold change\n(KCl - mannitol)", fontsize=BASE_FONT)
    ax.set_ylabel("Telescope subfamily-collapsed\nshrunken log2 fold change", fontsize=BASE_FONT)
    ax.set_title(
        "Six-hour effects across\nquantification approaches",
        fontsize=TITLE_FONT,
        fontweight="normal",
        pad=4,
    )
    ax.text(
        0.02,
        0.98,
        f"All: rho={rho:.2f}, P={p_value:.1e} (n=187)\n"
        f"IAP: rho={iap_rho:.2f}, P={iap_p_value:.1e} (n=19)",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=BASE_FONT,
        color=NAVY,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2},
        zorder=5,
    )
    ax.legend(
        frameon=False,
        loc="lower right",
        fontsize=BASE_FONT,
        borderaxespad=0.15,
        handletextpad=0.4,
        labelspacing=0.25,
    )
    style_axes(ax)
    panel_letter(fig, "E", 89, 86)
    return {
        "matched_targets": int(len(data)),
        "iap_related_targets": int(iap.sum()),
        "rho_all": float(rho),
        "p_all": float(p_value),
        "rho_iap": float(iap_rho),
        "p_iap": float(iap_p_value),
    }


def resolve_font_paths(font_dir: Path | None) -> tuple[Path, Path]:
    candidates = []
    if font_dir is not None:
        candidates.append(font_dir)
    candidates.extend([Path("/mnt/c/Windows/Fonts"), Path("C:/Windows/Fonts")])
    for directory in candidates:
        regular = directory / "arial.ttf"
        bold = directory / "arialbd.ttf"
        if regular.exists() and bold.exists():
            return regular, bold
    raise FileNotFoundError(
        "Arial font files were not found. Supply --font-dir containing arial.ttf and arialbd.ttf."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--font-dir", type=Path)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace a recognized prior output directory from this generator.",
    )
    args = parser.parse_args()

    global OUTDIR, PCA_SCORES, HEATMAP_VALUES, HEATMAP_STARS
    global COPY_PROPORTIONS, COPY_COMPOSITION, R34_TABLE, R42_TABLE, ANNOTATION_TABLE
    input_dir = args.input_dir.resolve()
    OUTDIR = args.output_dir.resolve()
    PCA_SCORES = input_dir / INPUT_FILENAMES["pca"]
    HEATMAP_VALUES = input_dir / INPUT_FILENAMES["heatmap_values"]
    HEATMAP_STARS = input_dir / INPUT_FILENAMES["heatmap_stars"]
    COPY_PROPORTIONS = input_dir / INPUT_FILENAMES["copy_proportions"]
    COPY_COMPOSITION = input_dir / INPUT_FILENAMES["copy_composition"]
    R34_TABLE = input_dir / INPUT_FILENAMES["dmr_null"]
    R42_TABLE = input_dir / INPUT_FILENAMES["telescope_de"]
    ANNOTATION_TABLE = input_dir / INPUT_FILENAMES["annotation"]

    inputs = [
        PCA_SCORES,
        HEATMAP_VALUES,
        HEATMAP_STARS,
        COPY_PROPORTIONS,
        COPY_COMPOSITION,
        R34_TABLE,
        R42_TABLE,
        ANNOTATION_TABLE,
    ]
    for path in inputs:
        if not path.exists():
            raise FileNotFoundError(path)
    if OUTDIR.exists() and any(OUTDIR.iterdir()):
        prior_summary = OUTDIR / "SUMMARY.json"
        existing_names = {path.name for path in OUTDIR.iterdir()}
        if existing_names <= {"inputs"}:
            pass
        elif not prior_summary.exists():
            raise SystemExit(f"Refusing to overwrite unrecognized output directory: {OUTDIR}")
        else:
            prior = json.loads(prior_summary.read_text())
            if prior.get("figure_id") != "FigS_TE_SalmonTE_family_timecourse" or not args.overwrite:
                raise SystemExit(
                    f"Recognized output exists at {OUTDIR}; pass --overwrite to replace it."
                )
    for subdir in ("figures",):
        (OUTDIR / subdir).mkdir(parents=True, exist_ok=True)

    arial_regular, arial_bold = resolve_font_paths(args.font_dir)
    for font_path in (arial_regular, arial_bold):
        fm.fontManager.addfont(str(font_path))
    matched_font = Path(fm.findfont("Arial", fallback_to_default=False))
    if matched_font.name.lower() not in {"arial.ttf", "arialbd.ttf"}:
        raise RuntimeError(f"Arial did not resolve correctly: {matched_font}")

    matplotlib.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": BASE_FONT,
            "axes.labelsize": BASE_FONT,
            "axes.titlesize": TITLE_FONT,
            "axes.titleweight": "normal",
            "xtick.labelsize": BASE_FONT,
            "ytick.labelsize": BASE_FONT,
            "legend.fontsize": BASE_FONT,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )

    fig = plt.figure(figsize=(PAGE_W_MM / 25.4, PAGE_H_MM / 25.4), facecolor="white")
    plot_panel_a(fig)
    plot_panel_b(fig)
    panel_c_summary = plot_panel_c(fig)
    panel_d_summary = plot_panel_d(fig)
    panel_e_summary = plot_panel_e(fig)

    pdf_path = OUTDIR / "figures/FigS_TE_SalmonTE_family_timecourse.pdf"
    png_path = OUTDIR / "figures/FigS_TE_SalmonTE_family_timecourse.png"
    fig.savefig(pdf_path, facecolor="white", metadata={"Title": "SalmonTE family-level supplementary figure"})
    fig.savefig(png_path, dpi=300, facecolor="white")
    plt.close(fig)

    print('Generated TE-family expression and genomic-context plots.')


if __name__ == '__main__':
    main()
