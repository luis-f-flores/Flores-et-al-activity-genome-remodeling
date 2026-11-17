"""Select the converged shared HMM with the highest balanced training likelihood."""
from pathlib import Path
import argparse
import json
import shutil
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--models', type=Path, required=True)
    parser.add_argument('--signals', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seeds', type=int, nargs='+', default=[12345, 23456, 34567])
    args = parser.parse_args()
    fits = [json.loads((args.models / f'seed_{seed}.metrics.json').read_text()) for seed in args.seeds]
    converged = [row for row in fits if row['converged']]
    if not converged:
        raise ValueError('No shared-HMM fit converged.')
    selected = max(converged, key=lambda row: row['final_balanced_loglik'])
    seed = selected['seed']
    args.output.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.models / f'seed_{seed}.model.json', args.output / 'selected_model.json')
    shutil.copyfile(args.models / f'seed_{seed}.decoded.npz', args.output / 'selected_decoded.npz')
    signals = dict(np.load(args.signals))
    states = np.load(args.output / 'selected_decoded.npz')['states']
    if states.shape != signals['condition_depth'].shape:
        raise ValueError('Decoded states and B1 measurements have different shapes.')
    np.savez_compressed(args.output / 'b1_signal_vectors.npz', **signals, states=states)
    (args.output / 'selected_fit.json').write_text(json.dumps(selected, indent=2) + '\n')
    print(f'Selected seed {seed} from {len(converged)} converged fits.')


if __name__ == '__main__':
    main()
