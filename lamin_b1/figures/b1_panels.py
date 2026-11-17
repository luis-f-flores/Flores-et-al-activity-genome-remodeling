"""Draw Lamin B1, RNA and IAP panels from the supplied measurement tables.

The figure selector covers Main 1–3 and S4–S8/S18; locus tracks use separate
renderers. Overlap annotations use
coverage-adjusted matching results, and GO panels use the supplied simplified
term lists and BH values."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
from scipy.cluster.hierarchy import dendrogram

COL = ["#0072B2", "#D55E00", "#4D4D4D"]
NAMES = ["TTX", "Bicuculline", "Bicuculline + nimodipine"]
SHORT = ["TTX", "Bic", "Bic + nimo"]
IDS = ["L9", "L12", "L11", "L14", "L16", "L10", "L13", "L15"]
GROUPS = [[0, 1], [2, 3, 4], [5, 6, 7]]
CONS = ["Bic_vs_TTX", "Bic_Nimo_vs_TTX", "Bic_Nimo_vs_Bic"]
HEAD = ["Bicuculline vs TTX", "Bicuculline + nimodipine vs TTX",
        "Bicuculline + nimodipine vs bicuculline"]
SAMPLES = ["Basal_A", "Basal_B", "Basal_C", "Bic_A", "Bic_B", "Bic_C",
           "Bic_nimo_A", "Bic_nimo_B", "Bic_nimo_C"]


class Inputs:
    def __init__(self, path, additional=()):
        self.path = Path(path)
        self.additional = [Path(p) for p in additional if p is not None]
        self.used = {}
        manifest = self.path / "inputs.json"
        self.hashes = json.loads(manifest.read_text())["files"] if manifest.exists() else {}
        self.overlap = None

    def file(self, name):
        path = next((p / name for p in [self.path, *self.additional] if (p / name).is_file()), self.path / name)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if name in self.hashes and digest != self.hashes[name]:
            raise ValueError(f"Input hash differs: {name}")
        self.used[name] = digest
        return path

    def table(self, name):
        return pd.read_csv(self.file(name), sep="\t")

    def json(self, name):
        if name == "DEG_LAD_matched_overlap.json" and self.overlap is not None:
            data = json.loads(self.overlap.read_text(encoding="utf-8"))
            return [dict(row, overlap=row.get("observed", row.get("overlap")))
                    for row in data["comparisons"]]
        return json.loads(self.file(name).read_text(encoding="utf-8"))

    def arrays(self, name):
        with np.load(self.file(name), allow_pickle=False) as data:
            return {k: data[k] for k in data.files}


def canvas(w, h):
    fig = plt.figure(figsize=(w / 25.4, h / 25.4))
    fig._mm = w, h
    return fig


def ax(fig, x, y, w, h):
    fw, fh = fig._mm
    a = fig.add_axes([x / fw, 1 - (y + h) / fh, w / fw, h / fh])
    a.tick_params(pad=2)
    return a


def text(fig, x, y, s, size=9, **kwargs):
    fw, fh = fig._mm
    return fig.text(x / fw, 1 - y / fh, s, fontsize=size,
                    va=kwargs.pop("va", "top"), **kwargs)


def title(fig, letter, label, x, y):
    text(fig, x, y, letter, 11, weight="bold")
    text(fig, x + 6, y + .3, label)


def legend(fig, y=3):
    fig.legend(handles=[Line2D([], [], color=c, lw=1.3, label=n)
                        for c, n in zip(COL, NAMES)], loc="upper center",
               bbox_to_anchor=(.5, 1 - y / fig._mm[1]), ncol=3,
               frameon=False, fontsize=8, columnspacing=.9, handlelength=1.2)


def save(fig, name, out):
    # Preserve component geometry while including labels extending past the axes.
    name = {"M01B": "Fig_1B", "M01E": "Fig_1E", "M01F": "Fig_1F", "M03AE": "Fig_3AE",
            "M02": "Fig_2", "S08ABC": "Fig_S08ABC"}.get(name, "Fig_"+name)
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight", pad_inches=.04,
                metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight", pad_inches=.04)
    plt.close(fig)


def venn(a, row):
    """Two circles with set-proportional areas and observed overlap area."""
    large = float(row["LAD"])
    small = float(row["DEG"])
    overlap = float(row["overlap"])
    r = math.sqrt(small / large)
    radius = 1.0

    def area(d):
        if d >= r + radius:
            return 0.0
        if d <= abs(radius - r):
            return math.pi * min(r, radius) ** 2
        return (r*r*math.acos((d*d+r*r-radius*radius)/(2*d*r))
                + radius*radius*math.acos((d*d+radius*radius-r*r)/(2*d*radius))
                - .5*math.sqrt((-d+r+radius)*(d+r-radius)*(d-r+radius)*(d+r+radius)))

    lo, hi = abs(radius-r), radius+r
    target = math.pi * overlap / large
    for _ in range(80):
        mid = (lo+hi)/2
        if area(mid) > target:
            lo = mid
        else:
            hi = mid
    d = (lo+hi)/2
    assert math.isclose(area(d)/math.pi*large, overlap, abs_tol=1e-7)
    a.add_patch(Circle((0, 0), radius, facecolor="#c4dce7", edgecolor="#7ca5b8", lw=.6, alpha=.78))
    a.add_patch(Circle((-d, 0), r, facecolor="#f5cdbb", edgecolor="#c7937b", lw=.6, alpha=.78))
    a.text(.15, 0, f"{int(large-overlap):,}", ha="center", va="center")
    a.text(.5, 1.13, "LAD genes", ha="center", fontsize=8)
    a.annotate(f"{int(small-overlap):,}", (-d-r*.55, 0), (-d-.5, -.85),
               ha="center", fontsize=8, arrowprops={"arrowstyle": "-", "lw": .5})
    a.annotate(f"{int(overlap):,}", ((-1-d+r)/2, 0), (-.9, 1.0),
               ha="center", fontsize=8, arrowprops={"arrowstyle": "-", "lw": .5})
    a.text(-d-.35, .65, "DEGs", ha="center", fontsize=8)
    statistic = (f"BH q = {row['BH_q']:.3g}" if row["BH_q"] is not None
                 else f"P = {row['P']:.2g}")
    a.text(0, -1.35, f"{row['fold']:.2f}-fold; matched {statistic}", ha="center", fontsize=8)
    a.set(aspect="equal", xlim=(min(-1.8, -d-r-.15), 1.3), ylim=(-1.5, 1.4))
    a.axis("off")


def go_panel(fig, rect, data, label, sizefactor=2.7):
    selected = data[data.retained_after_simplify].sort_values(["p.adjust", "ID"]).head(5)
    a = ax(fig, *rect)
    if selected.empty:
        a.axis("off")
        a.text(.5, .5, "No significant term", ha="center", transform=a.transAxes, fontsize=8)
        return
    ratios = np.array([float(s.split('/')[0])/float(s.split('/')[1]) for s in selected.GeneRatio])
    a.scatter(ratios, range(len(selected)), s=selected.Count*sizefactor,
              c=-np.log10(selected["p.adjust"]), cmap="Blues", norm=Normalize(1.3, 12),
              edgecolors="#43617a", lw=.35)
    a.set_yticks(range(len(selected)), selected.Description, fontsize=7)
    a.set(ylim=(len(selected)-.65, -.65), xlim=(0, max(.07, ratios.max()*1.4)))
    a.set_xlabel("Gene ratio", fontsize=8)
    a.set_title(label, fontsize=8)
    a.set_axisbelow(True)
    a.grid(axis="both", color="#dddddd", linewidth=.5)


def go_key(fig, x, y, sizefactor=2.7, counts=(5, 20)):
    from matplotlib.cm import ScalarMappable
    cb = fig.colorbar(ScalarMappable(norm=Normalize(1.3, 12), cmap="Blues"),
                     cax=ax(fig, x, y, 34, 2), orientation="horizontal", ticks=[2, 7, 12])
    cb.set_label("−log10(BH FDR)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    a = ax(fig, x+50, y-1, 44, 9)
    for xx, count in zip([.2, .65], counts):
        a.scatter(xx, .6, s=count*sizefactor, color="#5f99c6")
        a.text(xx+.09, .6, str(count), fontsize=8, va="center")
    a.set(xlim=(0, 1), ylim=(0, 1)); a.axis("off")
    text(fig, x+50, y+9, "Gene count", 8)


def m02(data, out):
    'Draw the transcriptome, 287-DEG, 173-gene, TPM and matched-overlap panels.'
    pca = data.table("RNA_PCA.tsv")
    vv = data.json("RNA_PCA_variance.json")
    hm = data.table("RNA_DEG_heatmap.tsv").sort_values("heatmap_row")
    assert len(hm) == 287 and hm.classification.tolist() == ["up"]*112 + ["down"]*175
    z = hm[[s+"_z" for s in SAMPLES]].to_numpy()
    np.testing.assert_allclose(z.mean(1), 0, atol=1e-13)
    trees = data.json("RNA_heatmap_linkage.json")
    for tree in trees:
        assert hm.loc[hm.classification == tree["direction"], "gene_id"].tolist() == tree["ordered_gene_ids"]
    scores = data.table("activity_scores_173.tsv")
    classes = ["rapidPRGs", "delayedPRGs", "SRGs"]
    amat = scores.pivot(index="gene_set", columns="sample", values="mean_gene_z_score").loc[classes, SAMPLES].to_numpy()
    assert amat.shape == (3, 9)
    fig = canvas(180, 150)
    legend(fig)
    title(fig, "A", "RNA PCA", 2, 10)
    a = ax(fig, 14, 23, 38, 27)
    for grp, col in zip(["TTX", "Bic", "Bic+Nimo"], COL):
        q = pca[pca.group == grp]
        a.scatter(q.pc1, q.pc2, s=21, color=col, edgecolors="white", linewidth=.4)
    a.set(xlim=(-33, 30), ylim=(-22, 25), xticks=[-30, 0, 30], yticks=[-20, 0, 20])
    a.set_xlabel(f"PC1 ({vv['PC1']:.1f}%)")
    a.set_ylabel(f"PC2 ({vv['PC2']:.1f}%)")
    title(fig, "B", "Bicuculline-responsive DEGs", 58, 10)
    a = ax(fig, 80, 23, 30, 32)
    im = a.imshow(z, aspect="auto", interpolation="none", cmap="RdBu_r", vmin=-2.5, vmax=2.5, rasterized=True)
    a.axhline(111.5, color="white", lw=.7)
    for xx in [2.5, 5.5]:
        a.axvline(xx, color="white", lw=.7)
    a.axis("off")
    bar = ax(fig, 80, 21, 30, 1.2)
    bar.imshow(np.repeat(np.arange(3), 3)[None, :], aspect="auto", cmap=ListedColormap(COL)); bar.axis("off")
    state = hm.TTX_LAD_display.map({"outside": 0, "unmapped": 1, "LAD": 2}).to_numpy()
    side = ax(fig, 77.5, 23, 1.1, 32)
    side.imshow(state[:, None], aspect="auto", cmap=ListedColormap(["white", "#bcbcbc", "#154968"]), vmin=0, vmax=2, interpolation="none"); side.axis("off")
    for direction, y0, hh in [("up", 23, 32*112/287), ("down", 23+32*112/287, 32*175/287)]:
        da = ax(fig, 72, y0, 4, hh)
        linkage = next(t["linkage"] for t in trees if t["direction"] == direction)
        dendrogram(np.asarray(linkage), orientation="left", no_labels=True, color_threshold=0, above_threshold_color="#777777", ax=da)
        da.invert_yaxis(); da.axis("off")
        for coll in da.collections: coll.set_linewidth(.4)
    text(fig, 70.5, 27, "112 up", 8, ha="right"); text(fig, 70.5, 45, "175 down", 8, ha="right")
    ea = ax(fig, 114, 23, 3.2, 32)
    ea.imshow(hm.Bic_Nimo_vs_Bic_log2FC.to_numpy()[:, None], aspect="auto", interpolation="none", cmap="RdBu_r", vmin=-2, vmax=2, rasterized=True)
    ea.axhline(111.5, color="white", lw=.7); ea.axis("off")
    text(fig, 119, 54, "Added nimodipine: log2FC", 7.5, rotation=90, va="bottom")
    for image, xx, ww, ticks in [(im, 80, 17, [-2.5, 0, 2.5]), (ea.images[0], 105, 12, [-2, 0, 2])]:
        cb = fig.colorbar(image, cax=ax(fig, xx, 61, ww, 1.4), orientation="horizontal", ticks=ticks)
        cb.ax.tick_params(labelsize=7, length=1, pad=1)
    text(fig, 80, 67, "Row z score", 8); text(fig, 106, 67, "log2FC", 8)
    title(fig, "C", "Activity response", 126, 10)
    a = ax(fig, 148, 25, 26, 25)
    imc = a.imshow(amat, aspect="auto", interpolation="none", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-1.25, vcenter=0, vmax=1.25))
    a.set_yticks(range(3), ["rPRGs (19)", "dPRGs (116)", "SRGs (38)"], fontsize=8)
    a.set_xticks([]); a.tick_params(length=0)
    for xx in [2.5, 5.5]: a.axvline(xx, color="white", lw=1)
    for spine in a.spines.values(): spine.set_visible(False)
    bar = ax(fig, 148, 22.5, 26, 1.2)
    bar.imshow(np.repeat(np.arange(3), 3)[None, :], aspect="auto", cmap=ListedColormap(COL)); bar.axis("off")
    cb = fig.colorbar(imc, cax=ax(fig, 150, 59, 22, 1.4), orientation="horizontal", ticks=[-1, 0, 1])
    cb.ax.tick_params(labelsize=7, length=1, pad=1); text(fig, 146, 65, "Mean gene z score", 8)
    title(fig, "D", "TTX RNA by B1 state", 2, 79)
    a = ax(fig, 22, 93, 48, 36)
    violins = data.json("TTX_TPM_state_violin.json")
    for i, (v, col) in enumerate(zip(violins, ["#abb3bc", "#78b9db", "#154968"]), 1):
        width, y = np.asarray(v["width"]), np.asarray(v["y"])
        a.fill_betweenx(y, i-width, i+width, color=col, alpha=.65)
        a.bxp([{"med": v["median"], "q1": v["q25"], "q3": v["q75"], "whislo": v["whislo"], "whishi": v["whishi"], "fliers": []}], positions=[i], widths=.12, showfliers=False, medianprops={"color": "black"}, boxprops={"lw": .5}, whiskerprops={"lw": .5}, capprops={"lw": .5})
    a.set(xticks=[1, 2, 3], xticklabels=["Non-B1", "Intermediate", "High"], ylim=(-.5, 13), yticks=[0, 6, 12], ylabel="TTX log2(TPM + 1)")
    a.tick_params(axis="x", rotation=35, labelsize=8)
    title(fig, "E", "DEG–LAD overlap", 91, 79)
    row = next(r for r in data.json("DEG_LAD_matched_overlap.json") if r["contrast"] == "Bic_vs_TTX")
    venn(ax(fig, 98, 91, 70, 50), row)
    save(fig, "M02", out)


def m01(data, out):
    fig = canvas(65, 65)
    title(fig, "B", "LAD-domain PCA", 2, 2)
    scores = data.table("domain_PCA_scores.tsv").set_index("sample")
    variance = data.table("domain_PCA_variance.tsv").percent
    a = ax(fig, 17, 14, 40, 35)
    for group, col in zip(GROUPS, COL):
        q = scores.loc[[IDS[j] for j in group]]
        a.scatter(q.PC1, q.PC2, c=col, s=18)
    a.set_xlabel(f"PC1 ({variance.iloc[0]:.1f}%)"); a.set_ylabel(f"PC2 ({variance.iloc[1]:.1f}%)")
    save(fig, "M01B", out)
    fig = canvas(180, 62)
    title(fig, "E", "Standardized B1 change in LAD segments", 2, 2)
    table = data.table("pairwise_B1_tests.tsv")
    palette = {"Not detected": "#BFC5C9", "Lower B1": "#0072B2", "Higher B1": "#D55E00"}
    for j, con in enumerate(CONS):
        a = ax(fig, 17+60*j, 19, 39, 28); d = table[table.contrast == con]
        for group in palette:
            q = d[d.B1_group == group]
            a.scatter(q.mean_z, q.B1_delta_z, c=palette[group], s=3 if group == "Not detected" else 5, alpha=.35 if group == "Not detected" else .8, lw=0, rasterized=True)
        a.axhline(0, color="#999999", ls="--", lw=.5)
        a.set(xlim=(-2, 4), ylim=(-5, 5), xticks=[-2, 0, 2, 4], yticks=[-4, 0, 4], xlabel="Mean B1 (z score)")
        a.text(.02, .97, f'{(d.B1_group=="Lower B1").sum()} lower / {(d.B1_group=="Higher B1").sum()} higher', transform=a.transAxes, va="top", fontsize=7)
        text(fig, 36+60*j, 10, HEAD[j].replace(" vs ", "\nvs "), 7.5, ha="center")
        if j == 0: a.set_ylabel("ΔB1 (z score)")
    save(fig, "M01E", out)
    fig = canvas(180, 66); title(fig, "F", "LAD-body, boundary and enrichment profiles", 2, 2)
    profiles = data.arrays("LAD_profiles.npz")
    for j, kind in enumerate(["body", "border"]):
        a = ax(fig, 17+60*j, 17, 40, 28)
        x = profiles[kind+"_x"]; mat = profiles[kind+"_means"][:, :8]
        for col, group in zip(COL, GROUPS):
            a.plot(x, mat[:, group], color=col, lw=.45, alpha=.25)
            a.plot(x, mat[:, group].mean(1), color=col, lw=1.2)
        a.set(ylim=(-.38, .32), yticks=[-.2, 0, .2]); a.axvline(0, color="#999", ls=":", lw=.5)
        if j == 0:
            a.axvline(1, color="#999", ls=":", lw=.5)
            a.set(xlim=(-.5, 1.5), xticks=[-.5, 0, 1, 1.5], xticklabels=["−500", "Start", "End", "+500"], xlabel="LAD body / flanks (kb)", ylabel="log2(B1/Input)")
        else: a.set(xlim=(-500, 500), xticks=[-500, 0, 500], xlabel="Boundary distance (kb)")
    a = ax(fig, 138, 17, 35, 28); ratio = data.table("calibration.tsv"); ratio = ratio[ratio.scheme == "condition_depth"]
    for j, name in enumerate(NAMES):
        y = ratio.loc[ratio.condition == name, "inside_outside_ratio"]
        a.scatter(j+np.linspace(-.08, .08, len(y)), y, s=18, color=COL[j]); a.plot([j-.2, j+.2], [y.median()]*2, color=COL[j], lw=1.2)
    a.axhline(1, color="#aaa", ls="--", lw=.5)
    a.set(xlim=(-.45, 2.45), ylim=(.98, 1.65), yticks=[1, 1.3, 1.6], xticks=[0, 1, 2], xticklabels=SHORT, ylabel="Inside / outside LADs")
    a.tick_params(axis="x", rotation=35, labelsize=8)
    save(fig, "M01F", out)


def m03(data, out):
    'Draw Figure 3A/E from the shared-HMM odds ratios and GO results.'
    fig = canvas(220, 151)
    title(fig, "A", "B1-change / DEG-direction associations", 2, 2)
    table = data.table("B1_RNA_OR.tsv")
    for j, con in enumerate(CONS[:2]):
        text(fig, 65+108*j, 12, HEAD[j].replace(" vs ", "\nvs "), 8, ha="center")
        a = ax(fig, 39+108*j, 27, 53, 21)
        for k, (bdir, rdir) in enumerate([("decreased", "up"), ("increased", "down")]):
            row = table[(table.contrast == con) & (table.scheme == "condition_depth") &
                        (table.B1_direction == bdir) & (table.DEG_direction == rdir)].iloc[0]
            y = 1-k
            a.plot([max(0, row.CI_low), min(4.5, row.CI_high)], [y]*2, color="#303030", lw=1)
            a.scatter([min(4.5, row.OR)], [y], s=18, color="black")
            if np.isfinite(row.CI_high) and row.CI_high > 4.5:
                a.plot(4.45, y, marker=">", color="#303030", ms=4)
            a.text(1.04, y, f"{row.OR:.2f}", transform=a.get_yaxis_transform(), va="center", fontsize=8)
            a.text(1.29, y, f"{row.BH_q:.2g}", transform=a.get_yaxis_transform(), va="center", fontsize=8)
        a.axvline(1, color="#999", lw=.6, ls="--")
        a.set(xlim=(0, 4.5), ylim=(-.5, 1.6), xticks=[0, 1, 2, 3, 4],
              yticks=[1, 0], yticklabels=["B1 ↓ / DEG ↑", "B1 ↑ / DEG ↓"], xlabel="Odds ratio")
        for xx, label in [(1.04, "OR"), (1.29, "q")]:
            a.text(xx, 1.4, label, transform=a.get_yaxis_transform(), fontsize=8)
        a.tick_params(labelsize=8)
    title(fig, "E", "GO cellular components: TTX-LAD DEGs", 2, 69)
    go = data.table("RNA_GO_CC.tsv")
    for j, direction in enumerate(["up", "down"]):
        rows = go[(go.contrast == "Bic_vs_TTX") & (go.direction == direction)]
        go_panel(fig, (57+108*j, 86, 35, 35), rows, direction.capitalize()+"regulated", 2.7)
    go_key(fig, 65, 136, 2.7, (5, 15))
    save(fig, "M03AE", out)


def three_sets(a, values, names, unit=""):
    centers = [(-.55, .24), (.55, .24), (0, -.62)]
    for xy, col in zip(centers, COL): a.add_patch(Circle(xy, .94, facecolor=col, edgecolor=col, alpha=.12, lw=.8))
    pos = {1: (-1.01, .4), 2: (1.01, .4), 4: (0, -1.18), 3: (0, .65), 5: (-.53, -.43), 6: (.53, -.43), 7: (0, -.05)}
    for k, xy in pos.items(): a.text(*xy, f"{values[k]:.1f}" if unit else f"{int(values[k]):,}", ha="center", va="center", fontsize=8)
    for x, y, name in [(-.8, 1.3, names[0]), (.8, 1.3, names[1]), (0, -1.8, names[2])]: a.text(x, y, name, ha="center", fontsize=8)
    a.set(aspect="equal", xlim=(-1.7, 1.7), ylim=(-2.05, 1.65)); a.axis("off")


def s04(data, out):
    fig = canvas(180, 168); legend(fig)
    pca = data.table("bin_PCA_scores.tsv"); pca = pca[(pca.scope == "all") & (pca.scheme == "condition_depth")].set_index("sample")
    variance = data.table("bin_PCA_variance.tsv"); variance = variance[(variance.scope == "all") & (variance.scheme == "condition_depth")].set_index("component").percent
    title(fig, "A", "Genome-wide B1 PCA", 2, 14); a = ax(fig, 22, 27, 51, 42)
    for group, col in zip(GROUPS, COL):
        q = pca.loc[[IDS[j] for j in group]]; a.scatter(q.PC1, q.PC2, s=24, color=col, edgecolors="white", linewidths=.35)
        for ii in group:
            r = pca.loc[IDS[ii]]; a.annotate(["A", "C", "A", "B", "C", "A", "B", "C"][ii], (r.PC1, r.PC2), xytext=(4, 4), textcoords="offset points", fontsize=8, color=col)
    a.margins(.19, .20); a.set_xlabel(f"PC1 ({variance.loc[1]:.1f}%)"); a.set_ylabel(f"PC2 ({variance.loc[2]:.1f}%)")
    corr = data.table("bin_correlations.tsv"); corr = corr[(corr.scheme == "condition_depth") & (corr.scope == "all")]
    matrix = corr.pivot(index="sample1", columns="sample2", values="Spearman").loc[IDS, IDS].to_numpy()
    calls = data.arrays("LAD_replication.npz"); states = calls["states"]; w = calls["clean_bp"]
    assert list(calls["IP_ids"]) == IDS and len(w) == 82945 and int(w.sum()) == 2448106300
    masks = states >= 1
    jac = np.array([[w[masks[:, i] & masks[:, j]].sum()/w[masks[:, i] | masks[:, j]].sum() for j in range(8)] for i in range(8)])
    for letter, label, rect, mat, low in [("B", "Sample correlation", (115, 27, 43, 43), matrix, .55), ("C", "LAD overlap", (20, 108, 43, 43), jac, min(.6, jac.min()))]:
        title(fig, letter, label, rect[0]-17, rect[1]-13)
        a = ax(fig, *rect); im = a.imshow(mat, vmin=low, vmax=1, cmap="Blues", interpolation="nearest")
        keys = ["T-A", "T-C", "B-A", "B-B", "B-C", "N-A", "N-B", "N-C"]
        a.set_xticks(range(8), keys, rotation=90, fontsize=7); a.set_yticks(range(8), keys, fontsize=7)
        for i in range(8):
            for j in range(8): a.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center", fontsize=6, color="white" if mat[i,j] > .84 else "#202020")
        fig.colorbar(im, cax=ax(fig, rect[0]+2, rect[1]+54, 38, 2), orientation="horizontal")
    condition = np.column_stack([masks[:, group].any(1) for group in GROUPS])
    sig = (condition*np.array([1, 2, 4])).sum(1)
    values = {k: float(w[sig == k].sum()/1e6) for k in range(8)}
    title(fig, "D", "Shared LAD territory (Mb)", 67, 95)
    three_sets(ax(fig, 70, 108, 48, 47), values, SHORT, "Mb")
    title(fig, "E", "LAD state composition", 124, 95)
    extent = data.table("condition_state_extent.tsv"); a = ax(fig, 140, 110, 35, 29)
    for j, name in enumerate(NAMES):
        q = extent[extent.condition == name].set_index("state"); iv = q.loc["Intermediate", "Mb"]; hv = q.loc["High", "Mb"]
        a.barh(2-j, iv, height=.56, color="#78b9db"); a.barh(2-j, hv, left=iv, height=.56, color="#154968")
    a.set(yticks=[2, 1, 0], yticklabels=SHORT, xlim=(0, 1650), xticks=[0, 750, 1500], xlabel="Territory (Mb)"); a.tick_params(labelsize=7)
    fig.legend(handles=[Patch(color="#78b9db", label="Intermediate"), Patch(color="#154968", label="High")], loc="upper center", bbox_to_anchor=(.83, .12), frameon=False, fontsize=8)
    save(fig, "S04", out)


def s05(data, out):
    fig = canvas(230, 235); legend(fig)
    calibration = data.table("calibration.tsv")
    for j, scheme in enumerate(["condition_depth", "condition_fly"]):
        title(fig, "AB"[j], ["Depth-normalized B1", "Spike-in calibrated B1"][j], 5+100*j, 14)
        a = ax(fig, 28+100*j, 27, 54, 36)
        for gi, name in enumerate(NAMES):
            q = calibration[(calibration.scheme == scheme) & (calibration.condition == name)]
            for (_, row), jitter in zip(q.iterrows(), np.linspace(-.075, .075, len(q))):
                xx = np.array([0, 1])+jitter; yy = [row.LAD_mean_log2, row.nonLAD_mean_log2]
                a.plot(xx, yy, color=COL[gi], lw=.65, alpha=.35); a.scatter(xx, yy, color=COL[gi], s=18, edgecolors="white", linewidths=.35)
        a.set(xlim=(-.3, 1.3), xticks=[0, 1], xticklabels=["LADs", "Non-LAD"], ylabel="Mean log2(B1/Input)")
        a.set_ylim((-.5, .4) if j == 0 else (-4, 3.5))
    title(fig, "C", "LAD genomic context", 5, 78)
    context = data.table("LAD_context.tsv"); classes = ["Shared LADs", "Variable LADs", "Outside all LADs"]
    for j, (key, label) in enumerate(zip(["scored_kb", "TSS_per_Mb", "AT_percent", "LINE_percent"], ["log10 length (kb)", "log10(TSS/Mb + 1)", "AT bases (%)", "LINE bases (%)"])):
        a = ax(fig, 40+47*j, 92, 37, 24)
        vals = [context.loc[context["class"] == cls, key].to_numpy() for cls in classes]
        vals = [np.log10(v) if j == 0 else np.log10(1+v) if j == 1 else v for v in vals]
        parts = a.violinplot(vals, positions=[2, 1, 0], vert=False, widths=.7, showextrema=False, points=200)
        for body, col in zip(parts["bodies"], ["#762A83", "#5E9FAD", "#ADADAD"]): body.set_facecolor(col); body.set_edgecolor(col); body.set_alpha(.65)
        for y, v in zip([2, 1, 0], vals):
            q = np.percentile(v, [25, 50, 75]); a.plot([q[0], q[2]], [y, y], color="white", lw=3); a.scatter(q[1], y, color="#222222", s=13, zorder=4)
        a.set(ylim=(-.6, 2.6), yticks=[2, 1, 0], yticklabels=classes if j == 0 else ["", "", ""], xlabel=label); a.tick_params(labelsize=7)
    title(fig, "D", "GO cellular components: differential B1", 5, 137)
    go = data.table("B1_GO_CC.tsv")
    for j, con in enumerate(CONS[:2]):
        for k, direction in enumerate(["Higher B1", "Lower B1"]):
            go_panel(fig, (73+112*j, 154+34*k, 35, 17), go[(go.contrast == con) & (go.direction == direction)], SHORT[j+1]+" vs TTX: "+direction)
    go_key(fig, 73, 222)
    save(fig, "S05", out)


def s06(data, out):
    d = data.table("RNA_diagnostic_genes.tsv"); full = data.json("RNA_DEG_counts.json")
    assert len(d) == 15954 and d.gene_id.is_unique
    x, y, z = d.BicvsTTX_log2FoldChange, d.Bic_NimovsTTX_log2FoldChange, d.Bic_NimovsBic_log2FoldChange
    direct = d.Bic_NimovsBic_DEG; extra = d.Bic_NimovsTTX_DEG & ~d.BicvsTTX_DEG
    miss_p = extra & (x.abs() >= .37) & (d.BicvsTTX_padj > .01)
    miss_both = extra & (x.abs() < .37) & (d.BicvsTTX_padj > .01)
    assert [int(q.sum()) for q in [direct, extra, miss_p, miss_both]] == [55, 5394, 3670, 1724]
    columns = [s+"_log2_normalized_count_plus1" for s in SAMPLES]
    for sl, key in [(slice(0, 3), "TTX_SD"), (slice(3, 6), "Bic_SD"), (slice(6, 9), "Bic+Nimo_SD")]:
        np.testing.assert_allclose(d[columns].iloc[:, sl].std(axis=1, ddof=1), d[key], rtol=1e-10, atol=1e-12)
    fig = canvas(180, 162); up, down = "#8656A4", "#3B8D6E"
    title(fig, "A", "Differentially expressed genes", 2, 3)
    a = ax(fig, 73, 16, 96, 28)
    for j, con in enumerate(["BicvsTTX", "Bic_NimovsTTX", "Bic_NimovsBic"]):
        r = full[con]; a.barh(2-j, r["down"], height=.5, color=down); a.barh(2-j, r["up"], left=r["down"], height=.5, color=up)
        a.text(r["DEGs"]+90, 2-j, f"{r['DEGs']:,}", va="center", fontsize=9)
    a.set(yticks=[2, 1, 0], yticklabels=["Bic / TTX", "Bic + nimo / TTX", "Bic + nimo / Bic"], xlim=(0, 6500), xticks=[0, 3000, 6000], xlabel="Number of DEGs")
    fig.legend(handles=[Patch(color=down, label="Downregulated"), Patch(color=up, label="Upregulated")], loc="upper right", frameon=False, fontsize=8, ncol=2)
    title(fig, "B", "Shared expression responses", 2, 57); a = ax(fig, 25, 72, 62, 62)
    a.hexbin(x, y, gridsize=55, mincnt=1, cmap="Greys", norm=LogNorm(), rasterized=True, linewidths=0)
    a.plot([-6.8, 6.8], [-6.8, 6.8], ls="--", lw=.7, c="#777777")
    for mask, col, label in [(direct & (z < 0), down, "48 downregulated"), (direct & (z > 0), up, "7 upregulated")]:
        a.scatter(x[mask], y[mask], s=9, color=col, edgecolor="white", lw=.2, zorder=3, label=label)
    a.set(xlim=(-6.8, 6.8), ylim=(-6.8, 6.8), xticks=[-6, 0, 6], yticks=[-6, 0, 6], xlabel="Bic / TTX log2 fold change", ylabel="Bic + nimo / TTX log2 fold change")
    a.text(.04, .97, f"r = {np.corrcoef(x,y)[0,1]:.3f}\nn = 15,954", transform=a.transAxes, va="top", fontsize=8)
    a.legend(frameon=False, fontsize=7, loc="lower right")
    title(fig, "C", "Bicuculline DEG criteria", 99, 57); text(fig, 106, 65, "Combined-treatment-only DEGs (n=5,394)", 7.5)
    a = ax(fig, 108, 79, 56, 20); counts = [int(miss_p.sum()), int(miss_both.sum())]
    a.barh([1, 0], counts, color=["#5E7E99", "#B2BEC8"], height=.66)
    for i, value in enumerate(counts): a.text(value+75, 1-i, f"{value:,}\n({value/5394:.1%})", va="center", fontsize=8)
    a.set(xlim=(0, 4900), yticks=[], xticks=[0, 2000, 4000], xlabel="Genes")
    for i, label in enumerate(["Effect passes;\nadjusted P misses", "Both criteria\nmissed"]): a.text(70, 1-i, label, va="center", fontsize=7, color="white" if i == 0 else "#202b33")
    title(fig, "D", "Replicate variability", 99, 116); a = ax(fig, 114, 129, 48, 20)
    for j, key in enumerate(["TTX_SD", "Bic_SD", "Bic+Nimo_SD"]):
        q = d[key].quantile([.25, .5, .75]).to_numpy(); a.errorbar(j, q[1], yerr=[[q[1]-q[0]], [q[2]-q[1]]], fmt="o", ms=4, lw=1.3, capsize=3, c=COL[j]); a.text(j, q[2]+.035, f"{q[1]:.3f}", ha="center", fontsize=8)
    a.set(ylim=(0, .7), xlim=(-.6, 2.6), yticks=[0, .3, .6], xticks=[0, 1, 2], xticklabels=SHORT, ylabel="Replicate SD"); a.tick_params(labelsize=8)
    save(fig, "S06", out)


def s07(data, out):
    fig = canvas(250, 224)
    title(fig, "A", "DEGs across RNA comparisons", 3, 4)
    rna = data.table("RNA_DEG_classification.tsv")
    sets = [set(rna.loc[(rna.contrast == c) & rna.classification.isin(["up", "down"]), "gene_id"]) for c in CONS]
    all_genes = set.union(*sets)
    values = {k: sum(sum((gene in s)*(1 << j) for j, s in enumerate(sets)) == k for gene in all_genes) for k in range(8)}
    three_sets(ax(fig, 7, 18, 62, 48), values, ["Bic/TTX", "+nimo/TTX", "+nimo/Bic"])
    rows = data.json("DEG_LAD_matched_overlap.json")
    for j, con in enumerate([CONS[2], CONS[1]]):
        title(fig, "BC"[j], HEAD[CONS.index(con)].replace("Bicuculline", "Bic"), 86+81*j, 4)
        venn(ax(fig, 88+81*j, 18, 70, 51), next(r for r in rows if r["contrast"] == con))
    go = data.table("RNA_GO_CC.tsv")
    for j, con in enumerate([CONS[2], CONS[1]]):
        title(fig, "DE"[j], "Directional GO cellular components", 3+125*j, 78)
        for k, direction in enumerate(["up", "down"]):
            go_panel(fig, (86+125*j, 92+39*k, 31, 22), go[(go.contrast == con) & (go.direction == direction)], direction.capitalize()+"regulated DEGs", .75)
    title(fig, "F", "Added nimodipine: B1 / DEG direction", 3, 177)
    table = data.table("B1_RNA_OR.tsv")
    q = table[(table.contrast == CONS[2]) & (table.scheme == "condition_depth")]
    a = ax(fig, 50, 192, 62, 18)
    for y, (bdir, rdir) in zip([1, 0], [("decreased", "up"), ("increased", "down")]):
        row = q[(q.B1_direction == bdir) & (q.DEG_direction == rdir)].iloc[0]
        if np.isfinite(row.OR):
            a.plot([max(0, row.CI_low), min(4.5, row.CI_high)], [y]*2, c="#303030", lw=1)
            a.scatter([min(4.5, row.OR)], [y], s=18, c="black")
            a.text(min(4.15, row.OR), y+.17, f"{row.OR:.2f}", ha="center", fontsize=8)
        else:
            a.text(3, y, "∞", fontsize=10, va="center")
        if np.isfinite(row.CI_high) and row.CI_high > 4.5:
            a.plot(4.45, y, marker=">", color="#303030", ms=4)
        a.text(1.05, y, f"{row.BH_q:.2g}", transform=a.get_yaxis_transform(), va="center", fontsize=8)
    a.axvline(1, color="#999", lw=.6, ls="--")
    a.set(xlim=(0, 4.5), ylim=(-.5, 1.55), xticks=[0, 1, 2, 3, 4], yticks=[1, 0], yticklabels=["B1 ↓ / DEG ↑", "B1 ↑ / DEG ↓"], xlabel="Odds ratio")
    a.text(1.05, 1.4, "q", transform=a.get_yaxis_transform(), fontsize=8)
    a.tick_params(labelsize=8)
    go_key(fig, 151, 195, .75, (50, 250))
    save(fig, "S07", out)


def s08(data, out):
    d = data.table("activity_B1_RNA.tsv"); rep = data.table("activity_B1_class_replicates.tsv")
    coverage = data.table("activity_LAD_coverage.tsv"); perm = data.table("activity_LAD_permutation.tsv")
    classes, sizes = ["rPRG", "dPRG", "SRG"], [19, 116, 38]
    assert len(d) == 173 and d.gene_name.is_unique and set(d.gene_name) == set(coverage.gene)
    blocks = [d[d.gene_set == cls].sort_values(["B1_TTX_mean", "gene_name"]) for cls in classes]
    assert [len(q) for q in blocks] == sizes
    order = pd.concat(blocks)
    for cls, size in zip(classes, sizes):
        q = coverage[coverage.gene_set == cls]; assert len(q) == size
        for key, name in zip(["TTX", "Bic", "Bic_Nimo"], NAMES):
            assert int(q[key+"_LAD75"].sum()) == int(perm[(perm.gene_set == cls) & (perm.condition == name)].LAD_genes.iloc[0])
        for ip in IDS:
            np.testing.assert_allclose(d.loc[d.gene_set == cls, "B1_"+ip].mean(), rep.loc[(rep.gene_set == cls) & (rep.IP == ip), "mean_body_log2"].iloc[0], atol=1e-12)
    fig = canvas(180, 189); title(fig, "A", "Activity-gene responses", 2, 2)
    cols_b = ["BicvsTTX_B1_delta", "Bic_NimovsTTX_B1_delta", "Bic_NimovsBic_B1_delta"]
    cols_r = ["BicvsTTX_log2FoldChange", "Bic_NimovsTTX_log2FoldChange", "Bic_NimovsBic_log2FoldChange"]
    for cols, x, width, cmap, lim, label in [(["B1_TTX_mean"], 24, 7, "Blues", (-1, 1), "TTX B1"), (cols_b, 39, 24, "PuOr_r", (-.5, .5), "ΔB1"), (cols_r, 72, 24, "RdBu_r", (-3, 3), "ΔRNA")]:
        text(fig, x+width/2, 13, label, 9, ha="center")
        cm = plt.get_cmap(cmap).copy(); cm.set_bad("#d9d9d9")
        a = ax(fig, x, 25, width, 105); im = a.imshow(np.ma.masked_invalid(order[cols].to_numpy(float)), aspect="auto", interpolation="nearest", cmap=cm, vmin=lim[0], vmax=lim[1]); a.set_yticks([]); a.set_xticks([])
        for cut in np.cumsum(sizes)[:-1]: a.axhline(cut-.5, color="white", lw=1.2)
        if len(cols) == 3:
            a.set_xticks(range(3), ["Bic − TTX", "Bic + nimo − TTX", "Bic + nimo − Bic"], rotation=65, ha="right", fontsize=8)
        cb = fig.colorbar(im, cax=ax(fig, x, 166, width, 2.3), orientation="horizontal", ticks=[lim[0], 0, lim[1]])
        cb.ax.tick_params(labelsize=7, pad=1, length=1); cb.set_label("log2(B1/Input)" if len(cols) == 1 else "Δlog2(B1/Input)" if cols == cols_b else "RNA log2FC", fontsize=7)
    cuts = np.r_[0, np.cumsum(sizes)]
    for i, (cls, size) in enumerate(zip(classes, sizes)): text(fig, 2, 25+105*(cuts[i]+size/2)/173-2, f"{cls}\n({size})", 8)
    title(fig, "B", "B1 response", 104, 2); a = ax(fig, 120, 34, 48, 46)
    for i, cls in enumerate(classes):
        for j, group in enumerate(GROUPS):
            y = rep.loc[(rep.gene_set == cls) & rep.IP.isin([IDS[k] for k in group]), "delta_from_TTX_mean"].to_numpy()
            xp = i+(j-1)*.22; a.scatter(xp+np.linspace(-.045, .045, len(y)), y, s=18, c=COL[j], lw=0); a.plot([xp-.09, xp+.09], [y.mean()]*2, c=COL[j], lw=1.1)
    a.axhline(0, color="#999", lw=.6, ls=":"); a.set(xlim=(-.5, 2.5), ylim=(-.2, .23), yticks=[-.2, 0, .2], xticks=[0, 1, 2], xticklabels=classes, ylabel="Δlog2(B1/Input)")
    for j, name in enumerate(SHORT): text(fig, 112, 14+5*j, name, 8, color=COL[j])
    title(fig, "C", "Entire-gene LAD coverage", 104, 95); a = ax(fig, 120, 112, 35, 49)
    for i, cls in enumerate(classes):
        for j, name in enumerate(NAMES):
            r = perm[(perm.gene_set == cls) & (perm.condition == name)].iloc[0]; yy = 8-(i*3+j)
            a.errorbar(r.null_mean, yy, xerr=[[r.null_mean-r.null_q025], [r.null_q975-r.null_mean]], fmt="|", c="#939ca2", elinewidth=1, capsize=2, markersize=4); a.scatter(r.observed_percent, yy, c=COL[j], s=20, zorder=3)
            text(fig, 167, 112+49*(8.8-yy)/9.6-1.5, f"{r.P_unadjusted:.3f}", 8, ha="center")
    text(fig, 167, 106, "P", 8, ha="center")
    a.set(xlim=(-3, 60), ylim=(-.8, 8.8), yticks=[7, 4, 1], yticklabels=classes, xticks=[0, 20, 40, 60], xlabel="Genes in LADs (%)")
    text(fig, 109, 176, "Gray: random-set central 95% interval", 7.5)
    save(fig, "S08ABC", out)


def s18(data, out):
    'Draw IAP components, including the median RNA change.'
    import csv
    def read(name):
        with data.file(name).open(encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))
    family=read("TE_family_effects.tsv")
    scatter=[r for r in family if r['new_measured']=='TRUE' and r['KCl_6h_measured']=='TRUE']
    x=np.array([float(r['KCl_Mann_6h']) for r in scatter]);y=np.array([float(r['Bic_TTX']) for r in scatter])
    def rank(v):
        _,i,n=np.unique(v,return_inverse=True,return_counts=True)
        return (np.cumsum(n)-(n-1)/2)[i]
    rho=np.corrcoef(rank(x),rank(y))[0,1]
    vals=read("IAP_display_values.tsv")
    dist=read("IAP_state_distribution.tsv");tests=read("IAP_enrichment_tests.tsv")
    CONDS=NAMES
    TARGETS=['IAPEY4_I','IAPA_MM','IAPEZI']
    IAPS=['IAPA_MM','IAPEZI','IAPEYI','IAPEY3_I','IAPEY4_I','IAPLTR1a_I_MM']
    fig=plt.figure(figsize=(180/25.4,235/25.4))
    def letter(s,x,y):fig.text(x,y,s,fontweight='bold',fontsize=11,va='top')
    # A: six-hour correlation, with labels for the families shown below.
    ax=fig.add_axes([.115,.715,.22,.215])
    mask=np.array([r['Name'] in IAPS for r in scatter])
    ax.axhline(0,color='#b7b7b7',lw=.6);ax.axvline(0,color='#b7b7b7',lw=.6)
    ax.scatter(x[~mask],y[~mask],s=9,color='#c8c8c8',linewidths=0)
    ax.scatter(x[mask],y[mask],s=19,color='#c46728',edgecolors='white',linewidths=.3,zorder=3)
    for name,offset in [('IAPA_MM',(-28,12)),('IAPEZI',(-20,15)),('IAPEY4_I',(7,7)),('IAPEYI',(8,-12)),('IAPEY3_I',(8,-3)),('IAPLTR1a_I_MM',(-8,-21))]:
     i=next(i for i,r in enumerate(scatter) if r['Name']==name)
     ax.annotate(name,(x[i],y[i]),xytext=offset,textcoords='offset points',fontsize=8,
      ha='right' if name=='IAPLTR1a_I_MM' else 'left',
      arrowprops=dict(arrowstyle='-',color='#777777',lw=.5),zorder=4)
    ax.set(xlim=(min(-1.55,x.min()-.1),max(3.1,x.max()+.1)),ylim=(min(-.87,y.min()-.1),max(.80,y.max()+.1)),
     xlabel='KCl − mannitol (6 h)\nlog2 mean-CPM ratio (+1)',ylabel='Bic − TTX\nlog2 mean-CPM ratio (+1)')
    ax.text(.04,.04,f'Spearman ρ = {rho:.3f}\nn = {len(scatter)} families',transform=ax.transAxes,fontsize=8,va='bottom')
    letter('A',.02,.977);fig.text(.07,.977,'TE-family responses',va='top',fontsize=9)
    # B: genomic composition of shared HMM states.
    colors=['#D8D8D8','#A78BBF','#51306E']
    ax=fig.add_axes([.435,.715,.245,.215])
    for gi,co in enumerate(CONDS):
     dr=sorted([r for r in dist if r['condition']==co],key=lambda r:int(r['state']))
     for yi,field in [(7-3*gi,'territory_bp'),(6-3*gi,'IAP_bp')]:
      v=np.array([float(r[field]) for r in dr]);v=100*v/v.sum();left=0
      for si,value in enumerate(v):
       ax.barh(yi,value,left=left,height=.65,color=colors[si],edgecolor='white',linewidth=.3)
       if value>=15:ax.text(left+value/2,yi,f'{value:.0f}',fontsize=8,ha='center',va='center',color='white' if si==2 else 'black')
       left+=value
     ax.text(0,7.65-3*gi,SHORT[gi],fontsize=8,va='bottom')
    ax.set(xlim=(0,100),ylim=(-.7,8.8),yticks=[7,6,4,3,1,0],yticklabels=['Genome','IAP']*3,
     xticks=[0,50,100],xlabel='Assessed sequence (%)')
    ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0,pad=3)
    letter('B',.35,.977);fig.text(.40,.977,'IAP distribution',va='top',fontsize=9)
    # C: spatial enrichment and six q values.
    q_values=[]
    for ei,(ep,title) in enumerate([('LAD / non-LAD','Inside / outside LADs'),('High / intermediate','High / intermediate')]):
     ax=fig.add_axes([.815,.83-.115*ei,.16,.075])
     for gi,co in enumerate(CONDS):
      row=next(r for r in tests if r['scope']=='all' and r['condition']==co and r['endpoint']==ep)
      v=float(row['density_ratio']);q_values.append(float(row['q_BH']))
      ax.plot(v,2-gi,'o',ms=3.5,color='#34383D')
      ax.text(v+.12,2-gi,f'{v:.2f}',va='center',fontsize=8)
     ax.axvline(1,color='#888888',ls='--',lw=.6)
     ax.set(xlim=(0,4.1),ylim=(-.5,2.5),yticks=[2,1,0],yticklabels=SHORT,xticks=[0,1,2,3])
     ax.tick_params(axis='y',length=0,pad=3);ax.spines['left'].set_visible(False)
     ax.set_title(title,fontsize=8,pad=6)
     if ei==0:ax.tick_params(axis='x',bottom=False,labelbottom=False);ax.spines['bottom'].set_visible(False)
     else:ax.set_xlabel('IAP density ratio')
    letter('C',.71,.977);fig.text(.76,.977,'IAP enrichment',fontsize=9,va='top')
    q_label=(f'All q = {q_values[0]:.4g}' if len(set(q_values))==1
             else f'q = {min(q_values):.3g}–{max(q_values):.3g}')
    fig.text(.845,.95,q_label,fontsize=8,ha='center',va='top')
    fig.legend(handles=[Patch(facecolor=c,label=s) for c,s in zip(colors,['Non-B1-associated','Intermediate','High'])],
     loc='center',bbox_to_anchor=(.64,.647),ncol=3,frameon=False,fontsize=8,handlelength=.9,columnspacing=.8)
    fig.legend(handles=[Line2D([],[],color='#c46728',marker='o',ls='',ms=4,label='Six highlighted IAP families')],
     loc='center',bbox_to_anchor=(.185,.647),frameon=False,fontsize=8,handletextpad=.3)
    # D–F: RNA; G–I: matched regional B1.
    plotted=[]
    for j,t in enumerate(TARGETS):
     left=.105+j*.315
     ra=fig.add_axes([left,.402,.245,.175]);ba=fig.add_axes([left,.105,.245,.19])
     mu=[];allv=[]
     for gi,co in enumerate(CONDS):
      rows=[r for r in vals if r['target']==t and r['assay']=='RNA' and r['condition']==co]
      assert len(rows)==3;plotted.extend(rows)
      v=np.array([float(r['value']) for r in rows]);mu.append(float(np.median(v)));allv.extend(v)
      ra.scatter(gi+np.linspace(-.065,.065,3),v,s=18,color='#34383D',edgecolors='white',linewidths=.35,zorder=3)
     ra.plot(range(3),mu,color='#34383D',lw=1.2);ym=max(allv)
     ra.set(xlim=(-.2,2.2),ylim=(0,ym*1.36),xticks=range(3),xticklabels=SHORT)
     ra.set_title(t,pad=8);letter(chr(68+j),left-.073,.607)
     if j==0:ra.set_ylabel('RNA (CPM)')
     for kind,label,col,ls,marker,shift in [('IAP-containing regions','IAP-containing regions','#8B3F8D','-','o',-.045),
      ('Matched non-IAP regions','Matched non-IAP regions','#60666D','--','s',.045)]:
      mu=[]
      for gi,co in enumerate(CONDS):
       rows=[r for r in vals if r['target']==t and r['assay']=='B1' and r['condition']==co and r['region_type']==kind]
       assert len(rows)==[2,3,3][gi];plotted.extend(rows)
       v=np.array([float(r['value']) for r in rows]);mu.append(v.mean())
       ba.scatter(gi+shift+np.linspace(-.05,.05,len(v)),v,s=17,color=col,marker=marker,edgecolors='white',linewidths=.3,zorder=3)
      ba.plot(range(3),mu,color=col,ls=ls,lw=1.5)
     ba.set(xlim=(-.2,2.2),ylim=(0,.56),xticks=range(3),xticklabels=SHORT,yticks=[0,.1,.2,.3,.4,.5])
     ba.set_title(t,pad=8);letter(chr(71+j),left-.073,.325)
     if j==0:ba.set_ylabel('B1 enrichment\nlog2(B1/Input)')
    assert sorted(json.dumps(r,sort_keys=True) for r in plotted)==sorted(json.dumps(r,sort_keys=True) for r in vals)
    fig.legend(handles=[Line2D([],[],color='#8B3F8D',marker='o',lw=1.5,ms=4,label='IAP-containing regions'),
     Line2D([],[],color='#60666D',marker='s',ls='--',lw=1.5,ms=4,label='Matched non-IAP regions')],
     loc='center',bbox_to_anchor=(.53,.045),ncol=2,frameon=False,fontsize=8.5)
    save(fig, "S18", out)


RENDERERS = {"M01": m01, "M02": m02, "M03": m03, "S04": s04, "S05": s05, "S06": s06, "S07": s07, "S08": s08, "S18": s18}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs')) / 'b1_current')
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--overlap", type=Path, help="Computed overlap RESULTS.json")
    parser.add_argument("--rna-inputs", type=Path, action="append", default=[], help="Generated RNA panel directory; may be repeated")
    parser.add_argument("--iap-inputs", type=Path, help="Generated IAP/B1 panel tables")
    parser.add_argument("--figure", choices=list(RENDERERS)+["all"], required=True)
    parser.add_argument("--font-dir", type=Path, help="Directory containing arial.ttf, arialbd.ttf and ariali.ttf; otherwise use Matplotlib's bundled DejaVu Sans")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("--output must be a new directory")
    plt.rcParams.update({"font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "pdf.fonttype": 42,
                         "axes.linewidth": .6, "axes.spines.top": False,
                         "axes.spines.right": False, "xtick.major.size": 2,
                         "ytick.major.size": 2})
    if args.font_dir:
        from matplotlib import font_manager
        for name in ["arial.ttf", "arialbd.ttf", "ariali.ttf"]:
            font_manager.fontManager.addfont(str(args.font_dir/name))
        plt.rcParams["font.family"] = "Arial"
    data = Inputs(args.inputs, [*args.rna_inputs, args.iap_inputs])
    data.overlap = args.overlap
    args.output.mkdir(parents=True)
    selected = list(RENDERERS) if args.figure == "all" else [args.figure]
    for key in selected:
        RENDERERS[key](data, args.output)
    receipt = {"figures": selected, "input_sha256": data.used,
               "outputs": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                           for p in sorted(args.output.iterdir()) if p.is_file()}}
    (args.output/"render_manifest.json").write_text(json.dumps(receipt, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"figures": selected, "output_files": list(receipt["outputs"])}))


if __name__ == "__main__":
    main()
