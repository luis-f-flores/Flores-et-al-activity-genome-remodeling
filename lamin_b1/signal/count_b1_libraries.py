"""Count the B1 IP and Input libraries with the autosomal and sex-chromosome masks."""
from pathlib import Path
import argparse
import subprocess
import sys

SAMPLES = ['L4', 'L6', 'L8', 'L2', 'L5', 'L7', 'L9', 'L12', 'L11', 'L14', 'L16', 'L10', 'L13', 'L15']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['bam-dir', 'autosomal-reference', 'sex-reference', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    for sample in SAMPLES:
        bam = args.bam_dir / f'{sample}.marked.bam'
        for name, reference, folder in [('count_unmasked_bins.py', args.autosomal_reference, 'autosomes'),
                                         ('count_sex_chromosome_fragments.py', args.sex_reference, 'sex_chromosomes')]:
            subprocess.run([sys.executable, str(Path(__file__).with_name(name)), 'count', '--sample', sample,
                            '--bam', str(bam), '--reference', str(reference), '--output', str(args.output / folder)], check=True)


if __name__ == '__main__':
    main()
