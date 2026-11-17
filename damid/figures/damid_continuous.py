#!/usr/bin/env python3
"""Calculate and render the numerical components of Fig. S01G-N.
"""
import os
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
os.environ['MPLBACKEND'] = 'Agg'

import argparse
import json
import importlib.util
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location('damid_continuous_analysis', ROOT / 'damid/continuous/measure_continuous_signal.py')
analysis = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(analysis)


def render(inputs, output, font_dir=None):
    inputs = Path(inputs)
    summary = json.loads((inputs / 'RESULTS.json').read_text())
    names = ['windows_2Mb.tsv.gz', 'window_read_support.tsv.gz', 'SD_points.tsv', 'rank_callouts.tsv', 'gene_windows.tsv']
    tables = {name: pd.read_csv(inputs / name, sep='\t', float_precision='round_trip') for name in names}
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    family = 'DejaVu Sans'
    if font_dir:
        for name in ('arial.ttf', 'arialbd.ttf', 'ariali.ttf'):
            font_manager.fontManager.addfont(str(Path(font_dir) / name))
        family = 'Arial'
    plt.rcParams.update({'font.family': family, 'font.size': 8, 'axes.labelsize': 8,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .65, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'xtick.major.size': 2.7, 'ytick.major.size': 2.7,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none', 'axes.unicode_minus': False})
    windows = tables['windows_2Mb.tsv.gz']
    windows = windows[windows.eligible].copy()
    support = tables['window_read_support.tsv.gz']
    points = tables['SD_points.tsv']
    labels = tables['rank_callouts.tsv']
    genes = tables['gene_windows.tsv']
    width, height = 170, 113
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))

    def axis(x, top, w, h):
        return fig.add_axes([x / width, 1 - (top + h) / height, w / width, h / height])

    def text(x, top, value, **kwargs):
        return fig.text(x / width, 1 - top / height, value, **kwargs)

    def letter(x, top, value):
        text(x, top, value, fontsize=11, fontweight='bold', va='center')

    blue, red = '#4475B6', '#D93027'
    heatmaps = {}
    for _, time in analysis.TIMES:
        ranked = windows[windows.time == time].sort_values('mean_change')
        order = pd.MultiIndex.from_frame(ranked[['chrom', 'start0', 'end0']])
        rows = []
        for field in ('fusion_count', 'Dam_count'):
            for condition, _, rep in analysis.PAIRS:
                sub = support[(support.time == time) & (support.condition == condition) & (support.replicate == rep)]
                values = sub.set_index(['chrom', 'start0', 'end0']).reindex(order)[field].to_numpy()
                if not np.isfinite(values).all() or values.min() < 10:
                    raise ValueError('Read heatmap is not aligned with the eligible ranks')
                rows.append(values)
        heatmaps[time] = np.vstack(rows)
    vmax = max(10000, 10 ** int(np.ceil(np.log10(max(a.max() for a in heatmaps.values())))))
    read_labels = ['Man A · B1', 'Man B · B1', 'KCl A · B1', 'KCl B · B1',
                   'Man A · Dam', 'Man B · Dam', 'KCl A · Dam', 'KCl B · Dam']
    for i, (_, time) in enumerate(analysis.TIMES):
        ax = axis(29 + i * 31, 11, 27, 22)
        image = ax.imshow(heatmaps[time], aspect='auto', interpolation='nearest', extent=(0, 100, 8, 0),
                          cmap='viridis', norm=LogNorm(vmin=10, vmax=vmax), rasterized=True)
        ax.set(yticks=np.arange(8) + .5, yticklabels=read_labels if i == 0 else [], xticks=[0, 50, 100])
        ax.tick_params(axis='y', length=0, pad=3, labelsize=6.7)
        ax.tick_params(axis='x', length=2, pad=2, labelsize=6.7)
        ax.axhline(4, color='white', lw=1.2)
        ax.set_title(time, fontsize=8, pad=3)
        ax.get_xticklabels()[0].set_ha('left')
        ax.get_xticklabels()[-1].set_ha('right')
        for spine in ax.spines.values():
            spine.set_visible(False)
    text(58, 40, 'Windows ranked by DamID change (%)', ha='center', fontsize=6.7)
    cb = fig.colorbar(image, cax=axis(89, 11, 1.4, 22), ticks=[10 ** p for p in range(1, int(np.log10(vmax)) + 1)])
    cb.ax.set_yticklabels([str(v) if v < 1000 else str(v // 1000) + 'k' for v in cb.get_ticks().astype(int)])
    cb.ax.tick_params(labelsize=6.5, length=2, pad=2)
    cb.ax.minorticks_off()
    text(90, 8, 'Reads', ha='center', fontsize=6.7)
    text(55, 3, 'Assigned reads per 2-Mb window', ha='center', fontsize=8)
    letter(2, 3, 'G')

    ax = axis(126, 11, 39, 22)
    positions = {('6 h', 'Mannitol'): 0, ('6 h', 'KCl'): .82,
                 ('Recovery', 'Mannitol'): 2.08, ('Recovery', 'KCl'): 2.90}
    for _, time in analysis.TIMES:
        means = {}
        for condition, color in (('Mannitol', blue), ('KCl', red)):
            sub = points[(points.time == time) & (points.condition == condition)].sort_values('replicate')
            values = sub.value.to_numpy()
            x, mean = positions[(time, condition)], float(values.mean())
            means[condition] = mean
            ax.errorbar(x, mean, yerr=values.std(ddof=1), fmt='none', ecolor=color, elinewidth=.8, capsize=2.3)
            ax.plot([x - .14, x + .14], [mean, mean], color=color, lw=1.2)
            ax.scatter(x + np.array([-.085, .085]), values, s=18, facecolor='white', edgecolor=color, lw=.8, zorder=4)
        center = (positions[(time, 'Mannitol')] + positions[(time, 'KCl')]) / 2
        ax.text(center, .87, f"{means['KCl'] / means['Mannitol']:.2f}", ha='center', va='center', fontsize=8)
        ax.text(center, -.12, time, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8)
    ax.set(xlim=(-.45, 3.35), ylim=(0, .92), xticks=[], yticks=[0, .2, .4, .6, .8],
           ylabel='SD of log$_2$\n(Dam–B1 / Dam)')
    ax.text(.5, 1.025, 'SD ratio (KCl / Mannitol)', ha='center', va='bottom', transform=ax.transAxes, fontsize=7)
    text(145, 3, 'Signal spread', ha='center', fontsize=8)
    letter(108, 3, 'H')

    for i, (_, time) in enumerate(analysis.TIMES):
        ax = axis(14 + i * 85, 50, 69, 24)
        ranked = windows[windows.time == time].sort_values('mean_change')
        ax.axvspan(0, 10, color='#EEEEEE', lw=0, zorder=0)
        ax.axvspan(90, 100, color='#EEEEEE', lw=0, zorder=0)
        ax.axhline(0, color='#BBBBBB', lw=.6, zorder=0)
        ax.plot(ranked.percentile, ranked.mean_change, color='#444444', lw=.85, zorder=2)
        for tail in ('Gain', 'Loss'):
            sub = labels[(labels.time == time) & (labels['tail'] == tail)].sort_values('percentile', ascending=False)
            for j, row in enumerate(sub.itertuples()):
                label = f'{row.chrom}:{row.start0 / 1e6:g}–{row.end0 / 1e6:g}  {row.example_genes}'
                position = (12, 2.72 - j * .82) if tail == 'Gain' else (12, -.72 - j * .82)
                ax.scatter(row.percentile, row.mean_change, s=11, color='#444444', zorder=4)
                ax.annotate(label, (row.percentile, row.mean_change), xytext=position, fontsize=6.7,
                            ha='left', va='center', bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': .6},
                            arrowprops={'arrowstyle': '-', 'color': '#555555', 'lw': .5,
                                        'relpos': (1, .5) if tail == 'Gain' else (0, .5)}, zorder=5)
        ax.set(xlim=(-1, 101), ylim=(-3.1, 3.2), xticks=[0, 10, 50, 90, 100], yticks=[-2, 0, 2], xlabel='Rank percentile')
        ax.set_title(time, fontsize=9, pad=4)
        if i == 0:
            ax.set_ylabel('Change in Dam–B1\nenrichment (log$_2$)', labelpad=3)
        letter(1 + i * 85, 46.5, 'IJ'[i])

    for j, gene in enumerate(('Bdnf', 'Dlg2', 'Cfap161', 'Stambpl1')):
        ax = axis(14 + j * 41.5, 92, 29.5, 15)
        ax.set_title(gene, fontstyle='italic', fontsize=9, pad=4)
        sub = genes[genes.gene == gene]
        for condition, color, offset in (('Mannitol', blue, -.055), ('KCl', red, .055)):
            means, sds = [], []
            for k, (_, time) in enumerate(analysis.TIMES):
                row = sub[sub.time == time].iloc[0]
                values = row[[condition + '_A', condition + '_B']].to_numpy(float)
                means.append(values.mean())
                sds.append(values.std(ddof=1))
                ax.scatter(k + offset + np.array([-.025, .025]), values, s=17,
                           facecolor='white', edgecolor=color, lw=.8, zorder=4)
            ax.errorbar(np.arange(2) + offset, means, yerr=sds, color=color, lw=1.05, elinewidth=.75,
                        capsize=2, marker='_', ms=8, mew=1.2, zorder=3)
        ax.axhline(0, color='#CCCCCC', lw=.6, zorder=0)
        ax.set(xticks=[0, 1], xticklabels=['6 h', 'Recovery'], xlim=(-.28, 1.28), ylim=(-1.7, 2.5), yticks=[-1, 0, 1, 2])
        if j == 0:
            ax.set_ylabel('Dam–B1 enrichment\n(log$_2$)', labelpad=3)
        letter(14 + j * 41.5, 88.5, 'KLMN'[j])
    fig.legend(handles=[Line2D([], [], color=blue, lw=1.2, label='Mannitol'),
                        Line2D([], [], color=red, lw=1.2, label='KCl')],
               loc='center', bbox_to_anchor=(.85, 1 - 40 / height), ncol=2, frameon=False,
               fontsize=8, handlelength=.8, columnspacing=.7, handletextpad=.4, borderpad=0)
    for ext in ('pdf', 'png'):
        fig.savefig(output / f'Fig_S01_G-N.{ext}', dpi=300, facecolor='white')
    plt.close(fig)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tables', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--font-dir', type=Path)
    args = parser.parse_args()
    render(args.tables, args.output, args.font_dir)


if __name__ == '__main__':
    main()
