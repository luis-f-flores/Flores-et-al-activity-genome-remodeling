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

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import time
import threading
from datetime import datetime, timezone
from pathlib import Path

TOOLS = {
    'bowtie2': 'bowtie2', 'build': 'bowtie2-build',
    'inspect': 'bowtie2-inspect', 'samtools': 'samtools',
    'trim': gb_path('@tool/trim_galore'),
    'cutadapt': gb_path('@tool/cutadapt'),
    'fastqc': gb_path('@tool/fastqc'),
    'curl': 'curl', 'multiqc': 'multiqc', 'pigz': 'pigz',
}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path, algorithm='sha256'):
    value = hashlib.new(algorithm)
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def save(path, value):
    temp = path.with_name(path.name + '.partial')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(path)


class Driver:
    def __init__(self, root, run):
        self.root, self.run = root.resolve(), run.resolve()
        self.refs = self.root / 'references'
        self.manifest_path = self.run / 'RUN_MANIFEST.json'
        self.manifest = json.loads(self.manifest_path.read_text())
        self.env = os.environ.copy()
        self.env['PYTHONNOUSERSITE'] = '1'
        self.env['PATH'] = os.environ.get('PATH', '')
        self.current_sample = None
        self.mutation_lock = threading.RLock()
        self.status('starting')

    def status(self, stage, **extra):
        save(self.run / 'STATUS.json', {'updated_utc': now(), 'pid': os.getpid(),
             'state': stage, 'sample': self.current_sample,
             'completed_samples': list(self.manifest.get('completed_samples', {})), **extra})
        print(now(), stage, self.current_sample or '', flush=True)

    def persist(self):
        self.manifest['updated_utc'] = now()
        save(self.manifest_path, self.manifest)

    def command(self, args, log, stdout_file=None):
        record = {'args': list(map(str, args)), 'started_utc': now(), 'log': str(log)}
        with self.mutation_lock:
            self.manifest.setdefault('commands', []).append(record)
            self.persist()
        with log.open('ab') as handle:
            if stdout_file:
                with stdout_file.open('wb') as output:
                    result = subprocess.run(record['args'], stdout=output, stderr=handle, env=self.env)
            else:
                result = subprocess.run(record['args'], stdout=handle, stderr=subprocess.STDOUT, env=self.env)
        with self.mutation_lock:
            record.update(finished_utc=now(), returncode=result.returncode)
            self.persist()
        if result.returncode:
            raise RuntimeError(f'Command failed ({result.returncode}); see {log}')

    def validate_raw(self):
        evidence = self.run / 'raw_integrity.json'
        if evidence.exists():
            record = json.loads(evidence.read_text())
            if all(Path(row['path']).stat().st_size == row['bytes'] and Path(row['path']).stat().st_mtime_ns == row['mtime_ns'] for row in record):
                return
            raise RuntimeError('Raw files changed since full integrity verification.')
        checks = []
        for row in self.manifest['raw_files']:
            self.current_sample = row['sample_id']
            self.status('raw_checksum', file=row['local_path'])
            path = Path(row['local_path'])
            before = path.stat()
            actual = sha(path, 'md5')
            after = path.stat()
            if before.st_size != row['expected_bytes'] or before.st_mtime_ns != after.st_mtime_ns or actual != row['expected_md5']:
                raise RuntimeError(f'Raw size/MD5 or stability check failed: {path}')
            self.status('raw_gzip_validation', file=str(path))
            self.command([TOOLS['pigz'], '-t', str(path)], self.run / 'logs/raw_integrity.log')
            checks.append({'sample_id': row['sample_id'], 'path': str(path), 'bytes': before.st_size,
                           'mtime_ns': before.st_mtime_ns, 'md5': actual, 'supplier_md5_matches': True,
                           'gzip_crc_passed': True})
        save(evidence, checks)
        self.manifest['raw_integrity'] = {'path': str(evidence), 'sha256': sha(evidence), 'files_verified': len(checks)}
        self.persist()
        self.current_sample = None

    def functional_check(self):
        self.status('tool_functional_check')
        with tempfile.TemporaryDirectory(prefix='tool_check_', dir=self.run / 'bam') as name:
            work = Path(name)
            seq = 'ACGT' * 25
            
            for mate in (1, 2):
                with gzip.open(work / f'test_{mate}.fq.gz', 'wt') as handle:
                    handle.write(f'@compatibility/{mate}\n{seq}\n+\n' + 'I' * 100 + '\n')
            self.command([TOOLS['trim'], '--paired', '--quality', '20', '--length', '30', '--cores', str(getattr(self, 'trim_cores', 2)),
                '--gzip', '--basename', 'compatibility', '--path_to_cutadapt', TOOLS['cutadapt'],
                '--output_dir', str(work), str(work / 'test_1.fq.gz'), str(work / 'test_2.fq.gz')], self.run / 'logs/tool_check.log')
            assert (work / 'compatibility_val_1.fq.gz').is_file()
            sam = work / 'test.sam'
            sam.write_text('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrTest\tLN:1000\n' +
                'test\t99\tchrTest\t101\t60\t50M\t=\t201\t150\t' + 'A' * 50 + '\t' + 'I' * 50 + '\n' +
                'test\t147\tchrTest\t201\t60\t50M\t=\t101\t-150\t' + 'T' * 50 + '\t' + 'I' * 50 + '\n')
            log = self.run / 'logs/tool_check.log'
            self.command([TOOLS['samtools'], 'sort', '-n', '-o', str(work / 'name.bam'), str(sam)], log)
            self.command([TOOLS['samtools'], 'fixmate', '-m', str(work / 'name.bam'), str(work / 'fix.bam')], log)
            self.command([TOOLS['samtools'], 'sort', '-o', str(work / 'coord.bam'), str(work / 'fix.bam')], log)
            self.command([TOOLS['samtools'], 'markdup', '-s', str(work / 'coord.bam'), str(work / 'marked.bam')], log)
            self.command([TOOLS['samtools'], 'index', str(work / 'marked.bam')], log)
            self.command([TOOLS['samtools'], 'quickcheck', '-v', str(work / 'marked.bam')], log)
            result = subprocess.run([TOOLS['samtools'], 'view', '-c', str(work / 'marked.bam')], capture_output=True, text=True, check=True)
            if result.stdout.strip() != '2':
                raise RuntimeError('SAMtools functional check did not retain both test records.')
        self.manifest['tool_functional_check'] = {'passed': True, 'utc': now(), 'log': str(self.run / 'logs/tool_check.log')}
        self.persist()

    def reference(self):
        self.refs.mkdir(exist_ok=True)
        manifest_path = self.refs / 'reference_manifest.json'
        prefix = self.refs / 'mm10_dm6'
        if manifest_path.exists():
            record = json.loads(manifest_path.read_text())
            if record.get('index_complete'):
                for path, item in record['retained_files'].items():
                    if not Path(path).is_file() or sha(path) != item['sha256']:
                        raise RuntimeError(f'Reference fingerprint differs: {path}')
                return prefix
        self.status('fly_reference_download')
        base = 'https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/'
        md5 = self.refs / 'dm6.ucsc_md5sum.txt'
        fly = self.refs / 'dm6.fa.gz'
        log = self.run / 'logs/reference_download.log'
        for filename, target in [('md5sum.txt', md5), ('dm6.fa.gz', fly)]:
            if not target.exists():
                partial = target.with_name(target.name + '.partial')
                self.command([TOOLS['curl'], '--fail', '--location', '--retry', '5', '--retry-delay', '5',
                    '--connect-timeout', '30', '--output', str(partial), base + filename], log)
                partial.replace(target)
        hashes = {line.split()[1].lstrip('*').removeprefix('./'): line.split()[0] for line in md5.read_text().splitlines() if line.strip()}
        if sha(fly, 'md5') != hashes['dm6.fa.gz']:
            raise RuntimeError('Drosophila FASTA does not match UCSC MD5.')
        mouse = Path(self.manifest['settings']['mouse_fasta'])
        self.status('combined_reference_preparation')
        combined = self.refs / 'mm10_dm6.fa'
        partial = self.refs / 'mm10_dm6.fa.partial'
        contigs, seen = [], set()
        with partial.open('wb') as output:
            for source, species, prefix_name in [(mouse, 'mm10', ''), (fly, 'dm6', 'dm6_')]:
                opener = gzip.open if source.suffix == '.gz' else open
                current, length = None, 0
                with opener(source, 'rb') as handle:
                    for line in handle:
                        if line.startswith(b'>'):
                            if current is not None:
                                contigs.append({'combined_name': prefix_name + current, 'source_name': current, 'species': species, 'length': length})
                            current = line[1:].split()[0].decode('ascii')
                            name = prefix_name + current
                            if name in seen:
                                raise RuntimeError(f'Duplicate reference contig: {name}')
                            seen.add(name)
                            length = 0
                            output.write(b'>' + name.encode('ascii') + b'\n')
                        else:
                            if current is None:
                                raise RuntimeError('FASTA sequence precedes first header.')
                            sequence = line.strip().upper()
                            length += len(sequence)
                            output.write(sequence + b'\n')
                    if current is not None:
                        contigs.append({'combined_name': prefix_name + current, 'source_name': current, 'species': species, 'length': length})
        expected = {line.split('\t')[0]: int(line.split('\t')[1]) for line in Path(self.manifest['settings']['mouse_fai']).read_text().splitlines()}
        actual = {row['source_name']: row['length'] for row in contigs if row['species'] == 'mm10'}
        if expected != actual:
            raise RuntimeError('Mouse source FASTA does not match the checked mm10 contig lengths.')
        partial.replace(combined)
        with (self.refs / 'contigs.tsv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['combined_name', 'source_name', 'species', 'length'], delimiter='\t')
            writer.writeheader()
            writer.writerows(contigs)
        record = {'created_utc': now(), 'mouse_source': str(mouse), 'mouse_sha256': sha(mouse),
            'fly_url': base + 'dm6.fa.gz', 'fly_build': 'UCSC dm6 / BDGP Release 6 + ISO1 MT, Aug 2014',
            'fly_ucsc_md5': hashes['dm6.fa.gz'], 'fly_md5_verified': True,
            'host_contigs': len(actual), 'fly_contigs': sum(r['species'] == 'dm6' for r in contigs),
            'contig_rule': 'Mouse names unchanged; all fly names prefixed dm6_. Sequence converted to uppercase, not hard-masked.',
            'index_prefix': str(prefix), 'index_complete': False}
        save(manifest_path, record)
        self.command([TOOLS['samtools'], 'faidx', str(combined)], self.run / 'logs/reference_faidx.log')
        self.status('combined_reference_index_build')
        self.command([TOOLS['build'], '--threads', '6', str(combined), str(prefix)], self.run / 'logs/reference_build.log')
        summary = self.refs / 'bowtie2_index_summary.txt'
        self.command([TOOLS['inspect'], '-s', str(prefix)], self.run / 'logs/reference_inspect.log', summary)
        observed = {line.split('\t')[1]: int(line.split('\t')[2]) for line in summary.read_text().splitlines() if line.startswith('Sequence-')}
        if observed != {row['combined_name']: row['length'] for row in contigs}:
            raise RuntimeError('Combined index contig names/lengths differ from combined FASTA.')
        indexes = list(self.refs.glob('mm10_dm6.*.bt2')) or list(self.refs.glob('mm10_dm6.*.bt2l'))
        if len(indexes) != 6:
            raise RuntimeError('Combined index does not contain six component files.')
        retained = [fly, md5, combined, Path(str(combined) + '.fai'), self.refs / 'contigs.tsv', summary] + indexes
        record['retained_files'] = {str(p): {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in retained}
        record['index_complete'] = True
        record['finished_utc'] = now()
        save(manifest_path, record)
        self.manifest['reference'] = {'path': str(manifest_path), 'sha256': sha(manifest_path), 'index_prefix': str(prefix)}
        self.persist()
        return prefix

    def fastqc(self, paths, sample, phase):
        destination = self.run / 'qc' / phase
        destination.mkdir(exist_ok=True)
        self.status('fastqc_' + phase)
        self.command([TOOLS['fastqc'], '--threads', str(getattr(self, 'qc_threads', 2)), '--outdir', str(destination)] + list(map(str, paths)),
                     self.run / 'logs' / f'{sample}.fastqc_{phase}.log')

    def sample(self, row, index):
        sample = row['sample_id']
        self.current_sample = sample
        completed = self.manifest.setdefault('completed_samples', {})
        if sample in completed:
            for path, record in completed[sample]['outputs'].items():
                if not Path(path).is_file() or Path(path).stat().st_size != record['bytes'] or sha(path) != record['sha256']:
                    raise RuntimeError(f'Completed checkpoint output changed: {path}')
            return
        r1, r2 = Path(row['r1_path']), Path(row['r2_path'])
        self.fastqc([r1, r2], sample, 'raw')
        self.status('trimming')
        trimmed = self.run / 'trimmed'
        self.command([TOOLS['trim'], '--paired', '--quality', '20', '--length', '30', '--cores', str(getattr(self, 'trim_cores', 2)), '--gzip',
            '--basename', sample, '--path_to_cutadapt', TOOLS['cutadapt'], '--output_dir', str(trimmed), str(r1), str(r2)],
            self.run / 'logs' / f'{sample}.trim.log')
        reads = [trimmed / f'{sample}_val_{mate}.fq.gz' for mate in (1, 2)]
        if not all(p.is_file() and p.stat().st_size for p in reads):
            raise RuntimeError(f'Missing paired trimmed output for {sample}')
        self.fastqc(reads, sample, 'trimmed')
        self.status('mapping_name_sort')
        bam = self.run / 'bam'
        name = bam / f'{sample}.name.pending.bam'
        coord = bam / f'{sample}.coord.pending.bam'
        pending = bam / f'{sample}.marked.pending.bam'
        final = bam / f'{sample}.marked.bam'
        align = [TOOLS['bowtie2'], '--very-sensitive', '--end-to-end', '--fr', '-I', '0', '-X', '1000', '--seed', '42',
                 '-p', str(getattr(self, 'alignment_threads', 6)), '--rg-id', sample, '--rg', f'SM:{sample}', '--rg', 'PL:ILLUMINA',
                 '-x', str(index), '-1', str(reads[0]), '-2', str(reads[1])]
        sort = [TOOLS['samtools'], 'sort', '-n', '-@', str(getattr(self, 'sort_threads', 2)), '-m', '1G', '-T', str(bam / f'{sample}.name_tmp'), '-o', str(name), '-']
        pipeline = {'commands': [align, sort], 'started_utc': now(), 'sample': sample}
        with self.mutation_lock:
            self.manifest.setdefault('pipelines', []).append(pipeline)
            self.persist()
        with (self.run / 'logs' / f'{sample}.bowtie2.log').open('ab') as alog, (self.run / 'logs' / f'{sample}.name_sort.log').open('ab') as slog:
            first = subprocess.Popen(align, stdout=subprocess.PIPE, stderr=alog, env=self.env)
            second = subprocess.Popen(sort, stdin=first.stdout, stderr=slog, env=self.env)
            first.stdout.close()
            rc_sort, rc_align = second.wait(), first.wait()
        with self.mutation_lock:
            pipeline.update(finished_utc=now(), returncodes=[rc_align, rc_sort])
            self.persist()
        if rc_align or rc_sort:
            raise RuntimeError(f'{sample}: alignment/name sorting failed ({rc_align}, {rc_sort}).')
        self.status('fixmate_coordinate_sort')
        fix = [TOOLS['samtools'], 'fixmate', '-m', '-u', str(name), '-']
        sort = [TOOLS['samtools'], 'sort', '-@', str(getattr(self, 'sort_threads', 2)), '-m', '1G', '-T', str(bam / f'{sample}.coord_tmp'), '-o', str(coord), '-']
        pipeline = {'commands': [fix, sort], 'started_utc': now(), 'sample': sample}
        with self.mutation_lock:
            self.manifest['pipelines'].append(pipeline)
            self.persist()
        with (self.run / 'logs' / f'{sample}.fixmate_sort.log').open('ab') as log:
            first = subprocess.Popen(fix, stdout=subprocess.PIPE, stderr=log, env=self.env)
            second = subprocess.Popen(sort, stdin=first.stdout, stderr=log, env=self.env)
            first.stdout.close()
            rc_sort, rc_fix = second.wait(), first.wait()
        with self.mutation_lock:
            pipeline.update(finished_utc=now(), returncodes=[rc_fix, rc_sort])
            self.persist()
        if rc_fix or rc_sort:
            raise RuntimeError(f'{sample}: fixmate/coordinate sorting failed ({rc_fix}, {rc_sort}).')
        self.status('duplicate_marking')
        self.command([TOOLS['samtools'], 'markdup', '-s', '-@', str(getattr(self, 'sort_threads', 2)), str(coord), str(pending)], self.run / 'logs' / f'{sample}.markdup.log')
        self.command([TOOLS['samtools'], 'quickcheck', '-v', str(pending)], self.run / 'logs' / f'{sample}.bam_validation.log')
        pending.replace(final)
        self.command([TOOLS['samtools'], 'index', '-@', str(getattr(self, 'sort_threads', 2)), str(final)], self.run / 'logs' / f'{sample}.bam_validation.log')
        self.status('bam_qc')
        outputs = [final, Path(str(final) + '.bai')] + reads
        for operation in ('flagstat', 'stats', 'idxstats'):
            output = self.run / 'qc' / f'{sample}.{operation}.txt'
            self.command([TOOLS['samtools'], operation, str(final)], self.run / 'logs' / f'{sample}.bam_validation.log', output)
            outputs.append(output)
        
        stat_text = (self.run / 'qc' / f'{sample}.stats.txt').read_text()
        count = re.search(r'^SN\traw total sequences:\t(\d+)', stat_text, re.M)
        align_text = (self.run / 'logs' / f'{sample}.bowtie2.log').read_text()
        pairs = re.findall(r'^(\d+) reads; of these:', align_text, re.M)
        if not count or not pairs or int(count.group(1)) != 2 * int(pairs[-1]):
            raise RuntimeError(f'{sample}: retained master BAM read count differs from Bowtie2 paired input count.')
        index_rows = list(csv.reader((self.run / 'qc' / f'{sample}.idxstats.txt').open(), delimiter='\t'))
        mouse_mapped = sum(int(r[2]) for r in index_rows if r[0] != '*' and not r[0].startswith('dm6_'))
        fly_mapped = sum(int(r[2]) for r in index_rows if r[0].startswith('dm6_'))
        checkpoint = {'completed_utc': now(), 'condition': row['condition'], 'assay': row['assay'],
            'trimmed_pairs_aligned': int(pairs[-1]), 'master_read_ends': int(count.group(1)),
            'mouse_mapped_read_ends_unfiltered': mouse_mapped, 'fly_mapped_read_ends_unfiltered': fly_mapped,
            'species_count_limit': 'Unfiltered mapped read ends; not eligible fragments or calibration factors.',
            'outputs': {str(p): {'bytes': p.stat().st_size, 'sha256': sha(p)} for p in outputs}}
        with self.mutation_lock:
            completed[sample] = checkpoint
            self.persist()
        
        for path in (name, coord):
            if path.parent.resolve() != bam.resolve():
                raise RuntimeError('Temporary cleanup path outside this run.')
            path.unlink()
        self.status('sample_complete')

    def execute(self):
        try:
            if shutil.disk_usage(self.run).free < 120 * 1024 ** 3:
                raise RuntimeError('Less than 120 GiB free at launch; check working storage before processing.')
            if not self.manifest.get('tool_functional_check', {}).get('passed'):
                self.functional_check()
            self.validate_raw()
            index = self.reference()
            rows = {row['sample_id']: row for row in self.manifest['sample_metadata']}
            for sample in self.manifest['queue']:
                self.sample(rows[sample], index)
            self.current_sample = None
            self.status('multiqc_summary')
            self.command([TOOLS['multiqc'], str(self.run / 'qc'), str(self.run / 'logs'), str(self.run / 'trimmed'),
                '--outdir', str(self.run / 'qc'), '--filename', 'multiqc_report', '--force'], self.run / 'logs/multiqc.log')
            self.manifest['finished_utc'] = now()
            self.manifest['state'] = 'complete_preprocessing'
            self.persist()
            self.status('complete_preprocessing', missing_expected_samples=['L1'],
                interpretation='Read trimming, alignment and quality checks complete.')
        except BaseException as exc:
            self.manifest['state'] = 'failed'
            self.manifest['error'] = repr(exc)
            self.persist()
            self.status('failed', error=repr(exc))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-root', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    Driver(args.package_root, args.run).execute()
