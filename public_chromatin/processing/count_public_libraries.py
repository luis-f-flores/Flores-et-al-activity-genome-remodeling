"""Count the deposited cortical ChIP and control libraries on a common mm9 grid."""
from pathlib import Path
import argparse
import csv
from count_public_chromatin_reads import count_reads


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--reads', type=Path, required=True)
    parser.add_argument('--sizes', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with args.samples.open(newline='', encoding='utf8') as handle:
        samples = list(csv.DictReader(handle, delimiter='\t'))
    if len({row['tag'] for row in samples}) != len(samples):
        raise ValueError('Library tags must be unique')
    sources = [(row['tag'], args.reads / row['filename']) for row in samples]
    for tag, source in sources:
        if not source.is_file():
            raise FileNotFoundError(source)
    for tag, source in sources:
        print('Counting', tag, flush=True)
        count_reads(source, args.sizes, args.output, tag)
