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
import collections
import contextlib
import csv
import fcntl
import gzip
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time
import zipfile

import pysam
from telescope_common import now, record, require, save, table, verify_records, version

PACKAGE = GB_WORK_ROOT / 'rna_package'
REFERENCE = PACKAGE / 'references/mm10_te_fragments'
TELESCOPE_PYTHON = sys.executable
SOURCE = Path(__file__).resolve().parent
CUTADAPT = gb_path('@tool/cutadapt')
ADAPTERS = ['AGATCGGAAGAGCACACGTCTGAACTCCAGTCA', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT']
IAP = ('IAPA_MM@', 'IAPEZI@', 'IAPEYI@', 'IAPEY3_I@')
MODULUS = 2**256
SPEC = dict(samples=['C6M', 'C6K'], k=500, threads=6, sort_threads=1,
            chunk_pairs=2000000, seed=20260912, reserve_gib=50,
            max_work_gib=100, target_rss_gib=24, max_rss_gib=28,
            min_available_gib=4, low_memory_seconds=10, swap_growth_gib=2,
            swap_pressure_seconds=60, ncpu=1, max_iter=100000,
            em_epsilon=1e-6, theta_prior=200000, pi_prior=0)


def module(name):
    spec = importlib.util.spec_from_file_location(name, SOURCE / (name + '.py'))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def environment():
    e = os.environ.copy()
    e.pop('PYTHONPATH', None)
    e.update(PYTHONNOUSERSITE='1', PYTHONHASHSEED='0', OMP_NUM_THREADS='1',
             OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
             TELESCOPE_ANNOTATION_CACHE=str(REFERENCE / 'telescope_annotation.pickle'))
    return e


@contextlib.contextmanager
def writer_lock(path):
    
    with Path(path).open('a+') as f:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        f.seek(0); f.truncate(); f.write(str(os.getpid())); f.flush()
        try:
            yield
        finally:
            f.seek(0); f.truncate(); f.flush()
            fcntl.flock(f, fcntl.LOCK_UN)


def configure_signals():
    def stop(signum, frame):
        raise InterruptedError(f'Signal {signum}: stop controller and its children')
    for signum in [signal.SIGTERM, signal.SIGINT, signal.SIGHUP]:
        signal.signal(signum, stop)


def assert_no_active_children(run):
    marker = Path(run)/'.active_children.json'
    groups = set()
    if marker.exists():
        groups = {int(x['pid']) for x in json.loads(marker.read_text()).get('children', [])}
    for p in Path('/proc').iterdir():
        if not p.name.isdigit() or int(p.name) == os.getpid():
            continue
        try:
            fields = (p/'stat').read_text().rsplit(')',1)[1].split()
            args = (p/'cmdline').read_bytes().split(b'\0')
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if fields[0] != 'Z' and int(fields[2]) in groups and any(a.startswith((str(run)+'/').encode()) for a in args):
            raise ValueError('A previous run child/process group is still active; inspect before resuming: '+p.name)


def memory():
    m = {s.split(':')[0]: int(s.split()[1]) for s in Path('/proc/meminfo').read_text().splitlines() if len(s.split()) >= 2}
    
    available = m.get('MemAvailable', m['MemFree']) / 1024**2
    return available, (m['SwapTotal'] - m['SwapFree']) / 1024**2


def tree_rss(pid):
    rows = []
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():
            continue
        try:
            fields = (p / 'stat').read_text().rsplit(')', 1)[1].split()
            rows.append((int(p.name), int(fields[1]), int(fields[21]) * os.sysconf('SC_PAGE_SIZE')))
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
    children = {pid}
    while True:
        more = {p for p, parent, rss in rows if parent in children}
        if more <= children:
            break
        children |= more
    
    return sum(rss for p, parent, rss in rows if p in children) / 2**30


class Guard:
    def __init__(self, run, spec, log):
        self.run, self.spec, self.log = Path(run), spec, Path(log)
        self.low_since = self.swap_since = None
        self.baseline_swap = memory()[1]
        self.last_check = self.last_disk = 0
        self.peak_rss = self.peak_work = self.peak_swap_growth = 0
        self.min_available = float('inf')

    def check(self, force=False):
        clock = time.monotonic()
        if not force and clock - self.last_check < 2:
            return
        self.last_check = clock
        available, swap = memory()
        rss = tree_rss(os.getpid())
        growth = max(0, swap - self.baseline_swap)
        self.peak_rss = max(self.peak_rss, rss)
        self.peak_swap_growth = max(self.peak_swap_growth, growth)
        self.min_available = min(self.min_available, available)
        self.low_since = (self.low_since or clock) if available < self.spec['min_available_gib'] else None
        self.swap_since = (self.swap_since or clock) if growth > self.spec['swap_growth_gib'] else None
        require(rss < self.spec['max_rss_gib'], 'Analysis RSS reached memory budget')
        require(available >= 2, 'Host available memory below emergency reserve')
        require(self.low_since is None or clock - self.low_since < self.spec['low_memory_seconds'], 'Sustained host memory pressure')
        require(self.swap_since is None or clock - self.swap_since < self.spec['swap_pressure_seconds'], 'Sustained growth in host swap usage')
        if force or clock - self.last_disk >= 10:
            self.last_disk = clock
            volumes = [self.run]
            if Path('/mnt/c').exists():
                volumes.append(Path('/mnt/c'))
            free = min(shutil.disk_usage(path).free for path in volumes) / 2**30
            used = sum(p.stat().st_size for p in self.run.rglob('*') if p.is_file() and not p.is_symlink()) / 2**30
            self.peak_work = max(self.peak_work, used)
            with self.log.open('a') as f:
                f.write(f'{now()}\t{rss:.4f}\t{available:.4f}\t{growth:.4f}\t{used:.4f}\t{free:.4f}\n')
            require(free >= self.spec['reserve_gib'], 'Working/backing volume free-space reserve reached')
            require(used < self.spec['max_work_gib'], 'New-run storage budget reached')

    def children(self, processes):
        save(self.run/'.active_children.json', dict(updated_utc=now(), controller_pid=os.getpid(),
             children=[dict(pid=p.pid,command=list(map(str,p.args))) for p in processes]))

    def summary(self):
        return dict(peak_analysis_rss_gib=self.peak_rss, min_host_available_gib=self.min_available,
                    peak_swap_growth_gib=self.peak_swap_growth, peak_run_gib=self.peak_work)


def pipeline(commands, logs, guard, allow_failure=False):
    processes, streams = [], []
    previous = None
    started = time.monotonic()
    try:
        for i, (command, logfile) in enumerate(zip(commands, logs)):
            f = Path(logfile).open('w'); streams.append(f)
            p = subprocess.Popen(list(map(str, command)), stdin=previous, stderr=f,
                                 stdout=subprocess.PIPE if i < len(commands)-1 else f,
                                 env=environment(), start_new_session=True)
            if previous is not None:
                previous.close()
            previous = p.stdout; processes.append(p)
            if hasattr(guard, "children"):
                guard.children(processes)
        while any(p.poll() is None for p in processes):
            guard.check()
            if any(p.poll() not in (None, 0) for p in processes):
                break
            time.sleep(1)
    finally:
        
        
        for p in processes:
            try:
                os.killpg(p.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        for p in processes:
            try:
                p.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(p.pid, signal.SIGKILL); p.wait()
            try:
                os.killpg(p.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        for f in streams:
            f.close()
        if hasattr(guard, "children"):
            guard.children([])
    codes = [p.returncode for p in processes]
    require(allow_failure or all(c == 0 for c in codes), f'Command failed ({codes}); inspect stage logs')
    guard.check(force=True)
    return dict(commands=commands, exit_codes=codes, elapsed_seconds=time.monotonic()-started,
                logs=[record(p) for p in logs])


def align_commands(first, second, output, threads=6):
    return [['bowtie2', '--very-sensitive-local', '-k', '500', '--score-min', 'L,0,1.6',
             '--seed', '20260912', '-p', str(threads), '-x', str(REFERENCE/'bowtie2/mm10'),
             '-1', str(first), '-2', str(second)],
            ['samtools', 'sort', '-n', '-@', '1', '-m', '512M', '-T', str(output)+'.sorttmp',
             '-o', str(output), '-']]


def load_command(bam, outdir, tag, ncpu=1):
    return [TELESCOPE_PYTHON, '-s', '-B', str(SOURCE/'telescope_disk_storage.py'), 'assign', str(bam), str(REFERENCE/'te_fragments.gtf'),
            '--attribute', 'locus', '--ncpu', str(ncpu), '--reassign_mode', 'exclude', '--conf_prob', '0.9',
            '--pi_prior', '0', '--theta_prior', '200000', '--em_epsilon', '1e-6', '--max_iter', '100000',
            '--skip_em', '--outdir', str(outdir), '--exp_tag', tag]


def fit_command(checkpoint, outdir, tag):
    return [TELESCOPE_PYTHON, '-s', '-B', '-m', 'telescope', 'resume', str(checkpoint), '--outdir', str(outdir),
            '--exp_tag', tag, '--max_iter', '100000', '--em_epsilon', '1e-6', '--theta_prior', '200000',
            '--pi_prior', '0', '--reassign_mode', 'exclude', '--conf_prob', '0.9']


def name_value(name):
    return int.from_bytes(hashlib.sha256(name).digest(), 'big')


def pairs(first, second):
    reader = module('fastq_pairs').fastq_record
    with gzip.open(first, 'rb') as a, gzip.open(second, 'rb') as b:
        while True:
            r1, r2 = reader(a, 1), reader(b, 2)
            if r1 is None and r2 is None:
                return
            require(r1 is not None and r2 is not None and r1[0] == r2[0], 'Missing/mismatched FASTQ mate')
            yield r1[0], r1[1], r2[1]


def split_pairs(first, second, folder, prefix, chunk_size, expected, guard):
    folder = Path(folder); chunks = []; total = 0
    input_hashes = [hashlib.sha256(), hashlib.sha256()]
    iterator = iter(pairs(first, second))
    while True:
        first_pair = next(iterator, None)
        if first_pair is None:
            break
        number = len(chunks) + 1
        paths = [folder/f'{prefix}.chunk{number:03d}_{m}.fq.gz' for m in (1,2)]
        n = digest = 0
        with contextlib.ExitStack() as stack:
            out = [stack.enter_context(gzip.GzipFile(fileobj=stack.enter_context(p.open('xb')),
                    mode='wb', filename='', mtime=0, compresslevel=1)) for p in paths]
            for name, one, two in itertools.chain([first_pair], itertools.islice(iterator, chunk_size-1)):
                out[0].write(one); out[1].write(two)
                input_hashes[0].update(one); input_hashes[1].update(two)
                digest = (digest + name_value(name)) % MODULUS
                n += 1; total += 1
                if total % 10000 == 0:
                    guard.check()
        chunks.append(dict(number=number, pairs=n, name_digest=f'{digest:064x}', reads=[record(p) for p in paths]))
    require(total == expected, f'Processed paired count mismatch: {total} versus {expected}')
    data = dict(total_pairs=total, chunk_size=chunk_size, chunks=chunks,
                uncompressed_sha256=[h.hexdigest() for h in input_hashes],
                membership='Count and additive SHA256 of names; merged BAM also requires exactly one primary record per mate per query name.')
    save(folder/f'{prefix}.chunks.json', data)
    return data


def bam_integrity(path, expected_pairs, expected_digest, guard, canonical=False):
    counts = collections.Counter(); digest = 0
    scientific = hashlib.sha256()
    with pysam.AlignmentFile(str(path), 'rb') as bam:
        require(bam.header.to_dict().get('HD', {}).get('SO') == 'queryname', 'Expected name-sorted BAM')
        reference = [line.split('\t') for line in (REFERENCE/'mm10.fa.fai').read_text().splitlines()]
        require(tuple(bam.references) == tuple(row[0] for row in reference) and
                tuple(bam.lengths) == tuple(int(row[1]) for row in reference),
                'BAM sequence dictionary differs from mm10.fa.fai')
        for name, group in itertools.groupby(bam, key=lambda r:r.query_name):
            primary = collections.Counter(); mapped = collections.Counter(); content = []
            for r in group:
                counts['alignment_records'] += 1
                require(r.is_paired and r.is_read1 != r.is_read2, 'Invalid paired alignment flags')
                mate = 1 if r.is_read1 else 2
                if not r.is_secondary and not r.is_supplementary:
                    primary[mate] += 1
                if not r.is_unmapped:
                    require(r.has_tag('AS'), 'Mapped alignment missing score')
                    mapped[mate] += 1
                if canonical:
                    v = r.to_string().split('\t')
                    content.append('\t'.join(v[:11]+sorted(v[11:])))
            require(primary == {1:1, 2:1}, f'Duplicate/missing primary mate at {name}')
            counts['pairs'] += 1
            counts['mapped_pairs'] += bool(mapped)
            counts['pairs_at_k500'] += max(mapped.values(), default=0) >= 500
            digest = (digest + name_value(name.encode())) % MODULUS
            if canonical:
                scientific.update(('\n'.join(sorted(content))+'\n').encode())
            if counts['pairs'] % 10000 == 0:
                guard.check()
    require(counts['pairs'] == expected_pairs, 'BAM lost or duplicated query-name groups')
    require(f'{digest:064x}' == expected_digest, 'BAM query names differ from processed reads')
    return dict(**counts, name_digest=f'{digest:064x}',
                canonical_alignment_sha256=scientific.hexdigest() if canonical else None)


def report(path):
    with Path(path).open() as f:
        line = next(f)
        require(line.startswith('## RunInfo'), 'Missing Telescope report header')
        info = dict(s.split(':',1) for s in line.rstrip().split('\t')[1:])
        rows = {}
        for row in csv.DictReader(f, delimiter='\t'):
            require(row['transcript'] not in rows, 'Duplicate Telescope locus')
            rows[row['transcript']] = row
    return info, rows


class Job:
    def __init__(self, run, manifest='RUN_MANIFEST.json'):
        self.run = Path(run); self.path = self.run/manifest
        self.state = json.loads(self.path.read_text())
        self.guard = Guard(run, self.state['spec'], self.run/'logs/resources.tsv')

    def write(self):
        self.state['updated_utc'] = now(); save(self.path, self.state)

    def stage(self, key, declared, work):
        previous = self.state['stages'].get(key)
        if previous and previous['status'] == 'complete':
            verify_records(previous['outputs'])
            return previous
        
        stamp = time.strftime('%Y%m%dT%H%M%SZ', time.gmtime())
        preserved = []
        files = list(map(Path, declared))
        for p in list(files):
            files.extend(p.parent.glob(p.name+'.sorttmp*'))
        for p in dict.fromkeys(files):
            require(p.resolve().is_relative_to(self.run.resolve()), 'Stage output escaped run directory')
            if p.exists():
                require(p.is_file() and not p.is_symlink(), 'Unexpected non-file stage output')
                dest = p.with_name(p.name+'.interrupted.'+stamp)
                require(not dest.exists(), 'Interrupted-output name already exists')
                p.rename(dest); preserved.append(record(dest))
        attempt = dict(status='running', started_utc=now(), preserved=preserved)
        if previous:
            attempt['previous_attempt'] = previous
        self.state['stages'][key] = attempt
        self.state.update(status='running', current_stage=key, pid=os.getpid()); self.write()
        started = time.monotonic()
        try:
            verify_records(self.state['code'])
            self.guard.check(force=True)
            result, outputs = work()
            self.guard.check(force=True)
            attempt.update(status='complete', result=result, outputs=[record(p) for p in outputs], completed_utc=now())
        except BaseException as e:
            attempt.update(status='failed', error=str(e), finished_utc=now())
            self.state.update(status='failed', error=str(e))
            raise
        finally:
            attempt['elapsed_seconds'] = time.monotonic()-started
            attempt['resources'] = self.guard.summary()
            self.write()
        print(now(), key, 'complete', flush=True)
        return attempt


def prepare(run, samples, reference, telescope_python, threads):
    require(not (run/'RUN_MANIFEST.json').exists(), 'Run already prepared')
    require(1 <= threads <= 32, 'Choose 1–32 alignment threads')
    run.mkdir(parents=True, exist_ok=True)
    (run/'logs').mkdir(exist_ok=True)
    with samples.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        required = {'sample','fastq1','fastq2','condition','replicate_label'}
        require(required <= set(reader.fieldnames or []), 'Sample table needs sample, fastq1, fastq2, condition and replicate_label columns')
        rows = list(reader)
    require(all(all(row.get(key) and row[key].strip() for key in required) for row in rows), 'Sample metadata fields must be nonempty')
    identifiers = [row['sample'] for row in rows]
    require(identifiers and len(identifiers) == len(set(identifiers)), 'Duplicate or missing sample IDs')
    require(all(re.fullmatch(r'[A-Za-z0-9_.-]+', sid) and sid not in {'.','..','logs','tables','source_reports'} for sid in identifiers), 'Invalid sample ID')
    inputs = {}
    for row in rows:
        sid = row['sample']
        (run/sid).mkdir(exist_ok=True)
        paths = [Path(row[f'fastq{mate}']).expanduser() for mate in (1,2)]
        paths = [p if p.is_absolute() else samples.parent/p for p in paths]
        inputs[sid] = [record(p.resolve()) for p in paths]
    reference = reference.resolve()
    refs = [reference/f for f in ['mm10.fa.fai','te_fragments.gtf','telescope_annotation.json','telescope_annotation.pickle']]
    index = sorted((reference/'bowtie2').glob('mm10*.bt2*'))
    require(len(index) == 6, 'Six Bowtie2 index files are required')
    sources = [record(SOURCE/f) for f in ['telescope_full_library.py','fastq_pairs.py',
               'telescope_annotation_cache.py','telescope_disk_storage.py','telescope_common.py']]
    state = dict(status='ready', created_utc=now(), spec=dict(SPEC,samples=identifiers,threads=threads),
        code=sources, inputs=inputs, sample_metadata=rows, references=[record(p) for p in refs+index], stages={},
        runtime=dict(reference=str(reference),telescope_python=os.path.abspath(os.path.expanduser(telescope_python))))
    save(run/'RUN_MANIFEST.json', state)


def execute(run):
    global REFERENCE, TELESCOPE_PYTHON
    runtime = json.loads((run/'RUN_MANIFEST.json').read_text())['runtime']
    REFERENCE = Path(runtime['reference'])
    TELESCOPE_PYTHON = runtime['telescope_python']
    with writer_lock(run/'.run.lock'):
        assert_no_active_children(run)
        job = Job(run); spec = job.state['spec']
        verify_records(job.state['code'] + job.state['references'])
        for items in job.state['inputs'].values():
            verify_records(items)
        for sid in spec['samples']:
            folder = run/sid
            clean = [folder/f'{sid}_{m}.fq.gz' for m in (1,2)]
            trim_json = folder/'preprocessing.json'; trim_log = run/'logs'/f'{sid}.preprocessing.log'
            def trim():
                cmd = [CUTADAPT, '-j','4','-a',ADAPTERS[0],'-A',ADAPTERS[1],'-e','0.1','-O','5',
                       '-q','20','--minimum-length','50','--max-n','0.1','--pair-filter','any','--json',str(trim_json),
                       '-o',str(clean[0]),'-p',str(clean[1])] + [x['path'] for x in job.state['inputs'][sid]]
                result = pipeline([cmd],[trim_log],job.guard)
                qc = json.loads(trim_json.read_text())['read_counts']
                require(qc['input'] > 0 and qc['output'] > 0, 'Preprocessing retained no read pairs')
                result['read_counts'] = qc
                return result, clean+[trim_json,trim_log]
            t = job.stage(sid+'.preprocess', clean+[trim_json,trim_log], trim)
            retained = t['result']['read_counts']['output']
            count = (retained+spec['chunk_pairs']-1)//spec['chunk_pairs']
            chunk_manifest = folder/f'{sid}.chunks.json'
            split_files = [folder/f'{sid}.chunk{c:03d}_{m}.fq.gz' for c in range(1,count+1) for m in (1,2)]
            def split():
                data = split_pairs(*clean,folder,sid,spec['chunk_pairs'],retained,job.guard)
                return data, split_files+[chunk_manifest]
            split_state = job.stage(sid+'.split', split_files+[chunk_manifest], split)
            chunks = split_state['result']['chunks']; bams = []
            for chunk in chunks:
                n = chunk['number']; bam = folder/f'{sid}.chunk{n:03d}.k500.name.bam'; bams.append(bam)
                logs = [run/'logs'/f'{sid}.chunk{n:03d}.{suffix}.log' for suffix in ['align','sort']]
                qcpath = folder/f'{sid}.chunk{n:03d}.bam_qc.json'
                def align(chunk=chunk,bam=bam,logs=logs,qcpath=qcpath):
                    result = pipeline(align_commands(*[x['path'] for x in chunk['reads']],bam,spec['threads']),logs,job.guard)
                    subprocess.run(['samtools','quickcheck',str(bam)],check=True)
                    qc = bam_integrity(bam,chunk['pairs'],chunk['name_digest'],job.guard)
                    save(qcpath,qc); result['bam_qc'] = qc
                    return result,[bam,qcpath]+logs
                job.stage(f'{sid}.align{n:03d}',[bam,qcpath]+logs,align)
            combined = folder/f'{sid}.k500.name.bam'; mq = folder/'combined_bam_qc.json'; ml = run/'logs'/f'{sid}.merge.log'
            def merge():
                cmd = ['samtools','merge','-n','-@','1','-c','-p',str(combined)]+list(map(str,bams))
                result = pipeline([cmd],[ml],job.guard)
                expected = sum(int(c['name_digest'],16) for c in chunks) % MODULUS
                qc = bam_integrity(combined,retained,f'{expected:064x}',job.guard)
                save(mq,qc); result['bam_qc'] = qc
                return result,[combined,mq,ml]
            job.stage(sid+'.merge',[combined,mq,ml],merge)
            tag = sid+'.k500'; checkpoint = folder/f'{tag}-checkpoint.npz'
            storage = folder/(tag+'.loading_storage.json'); loadlog = run/'logs'/f'{sid}.load.log'
            def load():
                result = pipeline([load_command(combined,folder,tag,spec['ncpu'])],[loadlog],job.guard)
                require(checkpoint.exists(),'No Telescope evidence checkpoint')
                with zipfile.ZipFile(checkpoint) as z:
                    require(z.testzip() is None,'Damaged checkpoint archive')
                return result,[checkpoint,storage,loadlog]
            spool = folder/f'{tag}-mappings.tmp'
            job.stage(sid+'.load',[checkpoint,storage,loadlog,spool],load)
            output = folder/f'{tag}-telescope_report.tsv'; fitlog = run/'logs'/f'{sid}.fit.log'
            def fit():
                result = pipeline([fit_command(checkpoint,folder,tag)],[fitlog],job.guard)
                match = re.search(r'EM (converged|terminated) after (\d+) iterations',fitlog.read_text())
                require(match and match[1]=='converged','Telescope did not converge; do not use counts')
                info, values = report(output)
                require(int(float(info['total_fragments']))==retained,'Telescope retained-pair total mismatch')
                result.update(iterations=int(match[2]),run_info=info,loci=len(values))
                return result,[output,fitlog]
            job.stage(sid+'.fit',[output,fitlog],fit)
        final(run,job)


def final(run,job):
    lines = []; summaries = {}; qc_rows = []
    reports = run/'source_reports'; reports.mkdir(exist_ok=True)
    (run/'tables').mkdir(exist_ok=True)
    for sid in job.state['spec']['samples']:
        info, values = report(run/sid/f'{sid}.k500-telescope_report.tsv')
        source_report = run/sid/f'{sid}.k500-telescope_report.tsv'
        with source_report.open('rb') as src, gzip.open(reports/f'primary_{sid}.telescope_report.tsv.gz','wb') as dst:
            shutil.copyfileobj(src,dst)
        metadata = next(row for row in job.state['sample_metadata'] if row['sample']==sid)
        biological = {key:row for key,row in values.items() if key!='__no_feature'}
        info.update(sample=sid,condition=metadata.get('condition',''),replicate_label=metadata.get('replicate_label',sid[0]),
            te_final_count_sum=sum(float(row['final_count']) for row in biological.values()),
            te_final_conf_sum=sum(float(row['final_conf']) for row in biological.values()))
        qc_rows.append(info)
        summaries[sid] = dict(fit=job.state['stages'][sid+'.fit']['result'],
                             bam=job.state['stages'][sid+'.merge']['result']['bam_qc'])
        for locus,row in values.items():
            if locus.startswith(IAP):
                lines.append([sid,locus,row['final_count'],row['final_conf'],row['unique_count']])
    table(run/'tables/sample_qc.tsv',list(qc_rows[0]),[[row.get(key,'') for key in qc_rows[0]] for row in qc_rows])
    sample_keys=list(job.state['sample_metadata'][0])
    table(run/'samples.tsv',sample_keys,[[row.get(key,'') for key in sample_keys] for row in job.state['sample_metadata']])
    table(run/'iap_counts.tsv',['sample','locus','assigned_fragments','confidence_ge_0.9','unique_assignment_fragments'],lines)
    save(run/'library_summary.json',dict(samples=summaries,alignment_alternatives=500))
    job.state.update(status='complete',current_stage=None,completed_utc=now());job.write()


def main():
    configure_signals()
    parser=argparse.ArgumentParser(description='Trim paired reads, align with Bowtie2 and quantify TE loci with Telescope.')
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--stage',choices=['prepare','run'],required=True)
    parser.add_argument('--samples',type=Path,help='TSV with sample, fastq1, fastq2, condition and replicate_label columns.')
    parser.add_argument('--reference-dir',type=Path,help='mm10.fa.fai, te_fragments.gtf, annotation cache and bowtie2/mm10 index.')
    parser.add_argument('--telescope-python',default=sys.executable)
    parser.add_argument('--threads',type=int,default=5)
    args=parser.parse_args(); run=args.run_dir.resolve()
    require(not run.is_relative_to(SOURCE.parents[1]),'Choose an output directory outside the code repository')
    if args.stage=='prepare':
        require(args.samples and args.reference_dir,'Preparation requires --samples and --reference-dir')
        prepare(run,args.samples.resolve(),args.reference_dir,args.telescope_python,args.threads)
    else:
        execute(run)


if __name__=='__main__':
    main()
