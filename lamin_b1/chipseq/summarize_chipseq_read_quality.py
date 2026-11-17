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

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
from datetime import datetime, timezone
import zipfile

ROOT = Path(gb_path('chipseq'))
RUN = ROOT / 'runs/2026-10-05_trim_map'
FASTQC = gb_path('@tool/fastqc')
MULTIQC = '/usr/bin/multiqc'


def now():
    return datetime.now(timezone.utc).isoformat()


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return {'path': str(path), 'bytes': path.stat().st_size,
            'sha256': digest.hexdigest()}


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def fastqc_name(path):
    name = path.name
    for suffix in ('.gz', '.fq', '.fastq'):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name + '_fastqc.zip'


def read_fastqc(path, phase):
    with zipfile.ZipFile(path) as archive:
        data_name = next(x for x in archive.namelist() if x.endswith('/fastqc_data.txt'))
        data = archive.read(data_name).decode()
        summary = archive.read(data_name.replace('fastqc_data.txt', 'summary.txt')).decode()
    filename = re.search(r'^Filename\t(.+)$', data, re.M).group(1)
    sid = re.match(r'^(L\d+)(?:_|\.)', filename).group(1)
    mate = int(re.search(r'_(1|2)\.(?:fq|fastq)(?:\.gz)?$', filename).group(1))
    total = int(re.search(r'^Total Sequences\t(\d+)$', data, re.M).group(1))
    modules = []
    for line in summary.splitlines():
        status, module, _ = line.split('\t')
        modules.append({'sample_id': sid, 'phase': phase, 'mate': mate,
                        'module': module, 'status': status, 'source': str(path)})
    return (sid, phase, mate), total, modules


def main():
    samples = list(csv.DictReader((ROOT / 'samples.tsv').open(), delimiter='\t'))
    samples = [row for row in samples if row['bam_path']]
    assert {r['sample_id'] for r in samples} == {f'L{i}' for i in range(2, 17)}
    qc = RUN / 'qc'
    if (qc / '14_QC_RESULT.json').exists():
        raise RuntimeError('Combined QC already exists; inspect its record before any repeat.')
    env = os.environ.copy()
    env['PYTHONNOUSERSITE'] = '1'
    env['PATH'] = str(Path(FASTQC).parent) + ':' + env['PATH']
    record = {'started_utc': now(), 'source': fingerprint(Path(__file__).resolve()),
              
              'commands': [], 'report_only': True}
    for tool in (FASTQC, MULTIQC):
        result = subprocess.run([tool, '--version'], env=env, capture_output=True, text=True, check=True)
        record.setdefault('software', {})[tool] = (result.stdout + result.stderr).strip()
    missing = []
    for row in samples:
        sid = row['sample_id']
        for phase, reads in [('raw', [Path(row['r1_path']), Path(row['r2_path'])]),
                             ('trimmed', [RUN / 'trimmed' / f'{sid}_val_{mate}.fq.gz' for mate in (1, 2)])]:
            absent = [path for path in reads if not (qc / phase / fastqc_name(path)).is_file()]
            if not absent:
                continue
            if not all(path.is_file() for path in absent):
                raise RuntimeError(f'Missing FASTQ needed for absent QC: {absent}')
            args = [FASTQC, '--threads', '2', '--outdir', str(qc / phase)] + [str(p) for p in absent]
            log = RUN / 'logs' / f'{sid}.fastqc_{phase}.log'
            command = {'args': args, 'log': str(log), 'started_utc': now()}
            print(f'Completing {sid} {phase} FastQC', flush=True)
            with log.open('a') as handle:
                subprocess.run(args, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
            command['completed_utc'] = now()
            record['commands'].append(command)
            missing.extend(str(qc / phase / fastqc_name(p)) for p in absent)
    counts, module_rows = {}, []
    evidence = {ROOT / 'samples.tsv'}
    for phase in ('raw', 'trimmed'):
        zips = sorted((qc / phase).glob('*.zip'))
        assert len(zips) == 30, (phase, len(zips))
        for path in zips:
            key, count, modules = read_fastqc(path, phase)
            assert key not in counts, key
            counts[key] = count
            module_rows.extend(modules)
            evidence.add(path)
    local = json.loads((RUN / 'RUN_MANIFEST.json').read_text())
    server_path = RUN / 'beast_transfer' / 'RUN_MANIFEST.json'
    server = json.loads(server_path.read_text())
    checkpoints = dict(local['completed_samples'])
    assert not set(checkpoints) & set(server['completed_samples'])
    checkpoints.update(server['completed_samples'])
    assert set(checkpoints) == {r['sample_id'] for r in samples}
    evidence.update([RUN / 'RUN_MANIFEST.json', server_path])
    rows = []
    for sample in samples:
        sid = sample['sample_id']
        stats_path = qc / f'{sid}.stats.txt'
        idx_path = qc / f'{sid}.idxstats.txt'
        bt_path = RUN / 'logs' / f'{sid}.bowtie2.log'
        evidence.update([stats_path, idx_path, bt_path, qc / f'{sid}.flagstat.txt',
                         RUN / 'logs' / f'{sid}.markdup.log', RUN / 'logs' / f'{sid}.trim.log'])
        stats = {}
        for line in stats_path.read_text().splitlines():
            if line.startswith('SN\t'):
                fields = line.split('\t')
                stats[fields[1].rstrip(':')] = float(fields[2])
        mouse = fly = 0
        for line in idx_path.read_text().splitlines():
            contig, _, mapped, _ = line.split('\t')
            if contig.startswith('dm6_'):
                fly += int(mapped)
            elif contig != '*':
                mouse += int(mapped)
        cp = checkpoints[sid]
        raw_pairs = counts[sid, 'raw', 1]
        trimmed_pairs = counts[sid, 'trimmed', 1]
        assert raw_pairs == counts[sid, 'raw', 2]
        assert trimmed_pairs == counts[sid, 'trimmed', 2]
        total = int(stats['raw total sequences'])
        mapped = int(stats['reads mapped'])
        assert total == 2 * trimmed_pairs == cp['master_read_ends']
        assert trimmed_pairs == cp['trimmed_pairs_aligned']
        assert mouse + fly == mapped
        assert mouse == cp['mouse_mapped_read_ends_unfiltered']
        assert fly == cp['fly_mapped_read_ends_unfiltered']
        bt = bt_path.read_text()
        assert int(re.findall(r'^(\d+) reads; of these:', bt, re.M)[-1]) == trimmed_pairs
        alignment = float(re.findall(r'([\d.]+)% overall alignment rate', bt)[-1])
        assert abs(alignment - 100 * mapped / total) <= 0.011
        inserts = []
        for line in stats_path.read_text().splitlines():
            if line.startswith('IS\t'):
                fields = line.split('\t')
                inserts.append((int(fields[1]), int(fields[2])))
        n = sum(count for _, count in inserts)
        def quantile(fraction):
            cumulative = 0
            for length, count in inserts:
                cumulative += count
                if cumulative >= fraction * n:
                    return length
        rows.append({'sample_id': sid, 'assay': sample['assay'], 'condition': sample['condition'],
                     'bio_rep': sample['bio_rep'], 'input_id': sample['input_id'],
                     'raw_pairs': raw_pairs, 'trimmed_pairs': trimmed_pairs,
                     'retained_pairs_pct': round(100 * trimmed_pairs / raw_pairs, 4),
                     'overall_alignment_pct': alignment, 'mouse_mapped_read_ends': mouse,
                     'fly_mapped_read_ends': fly, 'fly_pct_of_mapped_read_ends': round(100 * fly / mapped, 6),
                     'fly_pct_of_trimmed_read_ends': round(100 * fly / total, 6),
                     'duplicate_read_ends_pct': round(100 * stats['reads duplicated'] / total, 4),
                     'properly_paired_read_ends_pct': round(100 * stats['reads properly paired'] / total, 4),
                     'MQ0_mapped_read_ends_pct': round(100 * stats['reads MQ0'] / mapped, 4),
                     'insert_mean_bp': stats['insert size average'],
                     'insert_sd_bp': stats['insert size standard deviation'],
                     'insert_median_bp': quantile(.5), 'insert_p10_bp': quantile(.1), 'insert_p90_bp': quantile(.9)})
    table_path = qc / 'library_qc.tsv'
    modules_path = qc / 'fastqc_modules.tsv'
    write_tsv(table_path, rows)
    write_tsv(modules_path, module_rows)
    
    report_inputs = sorted(p for p in evidence if p.suffix == '.zip' or
                           (p.name.endswith(('.stats.txt', '.flagstat.txt', '.idxstats.txt', '.bowtie2.log', '.trim.log', '.markdup.log'))))
    file_list = qc / 'multiqc_inputs.txt'
    file_list.write_text('\n'.join(str(p) for p in report_inputs) + '\n')
    args = [MULTIQC, '--file-list', str(file_list), '--outdir', str(qc),
            '--filename', 'multiqc_report.html', '--title', 'Lamin B1 ChIP: 15 libraries, raw and trimmed QC',
            '--comment', 'Mouse mm10 + fly dm6; duplicate-retaining master BAMs. Fly counts measure aligned read ends.',
            '--module', 'fastqc', '--module', 'bowtie2', '--module', 'samtools', '--module', 'cutadapt',
            '--no-megaqc-upload', '--no-ansi']
    command = {'args': args, 'started_utc': now(), 'log': str(RUN / 'logs/14_combined_qc.log')}
    print('Generating combined MultiQC 1.12', flush=True)
    with (RUN / 'logs/14_combined_qc.log').open('a') as handle:
        subprocess.run(args, env=env, stdout=handle, stderr=subprocess.STDOUT, check=True)
    command['completed_utc'] = now()
    record['commands'].append(command)
    assert (qc / 'multiqc_report.html').stat().st_size > 1000
    record.update({'completed_utc': now(), 'state': 'complete', 'libraries': 15,
                   'raw_fastqc_mates': 30, 'trimmed_fastqc_mates': 30,
                   'completed_missing_fastqc': missing, 'count_reconciliation_passed': True,
                   'fly_unit': 'Mapped read ends; includes duplicates, singletons and all mapping qualities.',
                   'fly_percent_denominator': 'All mapped mouse + fly read ends.',
                   'paired_fragment_calibration': 'This report summarizes read-end QC; fragment-based calibration is calculated separately.',
                   'input_files': [fingerprint(p) for p in sorted(evidence)],
                   'outputs': [fingerprint(p) for p in [table_path, modules_path, file_list, qc / 'multiqc_report.html']]})
    (qc / '14_QC_RESULT.json').write_text(json.dumps(record, indent=2) + '\n')
    handoff_path = RUN / 'beast_handoff.json'
    handoff = json.loads(handoff_path.read_text())
    handoff['combined_local_qc'] = {'completed_utc': record['completed_utc'],
        'state': 'complete', 'source': str(Path(__file__).resolve()),
        'result': str(qc / '14_QC_RESULT.json'), 'report': str(qc / 'multiqc_report.html'),
        'table': str(table_path), 'calibration_performed': False, 'sample_exclusions': []}
    temporary = handoff_path.with_suffix('.qc_pending')
    temporary.write_text(json.dumps(handoff, indent=2) + '\n')
    temporary.replace(handoff_path)
    print(json.dumps(rows, indent=2), flush=True)
    print('COMPLETE: 15-library QC; 60 FastQC mates; counts reconcile with source checkpoints.', flush=True)


if __name__ == '__main__':
    main()
