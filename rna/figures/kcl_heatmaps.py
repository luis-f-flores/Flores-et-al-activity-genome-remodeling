"""Calculate and draw the gene and activity-class TPM heatmaps."""
from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'expression'))
from kcl_tpm import calculate, COLS, CLASSES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    z, activity = calculate(args.inputs, args.output/'tables')
    order = pd.read_csv(args.inputs/'TPM_gene_order.csv')
    plt.rcParams.update({'font.size': 9, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 4.2), gridspec_kw={'width_ratios': [1, 1]}, layout='constrained')
    for ax, values, limits, title in zip(axes, (z, activity), ((-2.5, 2.5), (-2, 2)), ('Differentially expressed genes', 'Activity-response classes')):
        im = ax.imshow(values, aspect='auto', interpolation='nearest', cmap='RdBu_r', vmin=limits[0], vmax=limits[1], rasterized=True)
        ax.set_title(title)
        ax.set_xticks(range(len(COLS)), COLS, rotation=90, fontsize=7)
        fig.colorbar(im, ax=ax, shrink=.6, label='Row Z score of log2(TPM + 1)')
    axes[0].set_yticks([])
    boundaries = np.flatnonzero(order.Cluster.to_numpy()[1:] != order.Cluster.to_numpy()[:-1])+.5
    for y in boundaries: axes[0].axhline(y, color='black', linewidth=.4)
    axes[1].set_yticks(range(len(CLASSES)), CLASSES)
    fig.savefig(args.output/'Fig_5.pdf')
    fig.savefig(args.output/'Fig_5.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
