"""Count nonmitochondrial Drosophila fragments used for B1 spike-in scaling."""
from pathlib import Path
import argparse
import csv
import pysam

IPS = ['L9', 'L12', 'L11', 'L14', 'L16', 'L10', 'L13', 'L15']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bam-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for sample in IPS:
        count = 0
        with pysam.AlignmentFile(str(args.bam_dir / f'{sample}.marked.bam'), 'rb', threads=1) as bam:
            for contig in bam.references:
                if not contig.startswith('dm6_') or contig == 'dm6_chrM':
                    continue
                for read in bam.fetch(contig):
                    if (read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate
                            or read.is_qcfail or not read.is_read1 or not read.is_paired or not read.is_proper_pair
                            or read.mate_is_unmapped or read.next_reference_id != read.reference_id
                            or not read.template_length or not read.has_tag('MQ')):
                        continue
                    if min(read.mapping_quality, read.get_tag('MQ')) >= 10:
                        count += 1
        rows.append(dict(sample_id=sample, fly_nonmitochondrial_nonduplicate_proper_pairs_MAPQ10=count))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
