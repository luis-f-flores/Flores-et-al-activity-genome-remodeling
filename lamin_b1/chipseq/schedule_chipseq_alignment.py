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
import concurrent.futures
import fcntl
import hashlib
import importlib.util
import json
import os
import signal
import shutil
import subprocess
import threading
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path

GIB = 1024 ** 3


def now():
    return datetime.now(timezone.utc).isoformat()


def save(path, value):
    temporary = path.with_name(path.name + '.partial')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def processes():
    rows = {}
    for path in Path('/proc').glob('[0-9]*/stat'):
        try:
            text = path.read_text()
            fields = text[text.rfind(')') + 2:].split()
            rows[int(path.parent.name)] = {'state': fields[0], 'ppid': int(fields[1]),
                'rss': int(fields[21]) * os.sysconf('SC_PAGE_SIZE'),
                'exit_code': int(fields[49]) if len(fields) > 49 else None}
        except (OSError, ValueError):
            continue
    return rows


def descendants(pid, rows):
    children = set()
    frontier = {pid}
    while frontier:
        following = {n for n, row in rows.items() if row['ppid'] in frontier} - children
        children.update(following)
        frontier = following
    return children


def available_memory():
    fields = {line.split(':', 1)[0]: int(line.split(':', 1)[1].split()[0]) * 1024
              for line in Path('/proc/meminfo').read_text().splitlines()}
    return fields.get('MemAvailable', fields['MemFree'])


def handoff_ready(state, completed, first, wait_library):
    boundary = state['state'] in {'sample_complete', 'fastqc_raw'} or state.get('sample') in completed
    return first in completed and wait_library in completed and boundary and state['state'] != 'failed'


class Coordinator:
    def __init__(self, root, run):
        self.root, self.run = root, run
        self.plan_path = run / 'parallel_plan.json'
        self.plan = json.loads(self.plan_path.read_text())
        self.manifest_path = run / 'RUN_MANIFEST.json'
        self.lock = threading.RLock()
        self.worker_stages = {}
        self.manifest = None

    def handoff(self):
        
        receipt_path = self.run / 'parallel_takeover_QC.json'
        receipt = json.loads(receipt_path.read_text())
        if not receipt.get('passed'):
            raise RuntimeError('Verified stop/checkpoint receipt is required.')
        rows = processes()
        for pid in receipt['old_process_ids']:
            if pid in rows and rows[pid]['state'] not in {'Z','X'}:
                raise RuntimeError('Old serial process is still present: ' + str(pid))
        for path, fingerprint in receipt['verified_files'].items():
            stat = Path(path).stat()
            if stat.st_size != fingerprint['bytes'] or stat.st_mtime_ns != fingerprint['mtime_ns']:
                raise RuntimeError('Verified checkpoint file changed: ' + path)
        self.manifest = json.loads(self.manifest_path.read_text())
        for sample, original in receipt['completed_records'].items():
            if self.manifest['completed_samples'].get(sample) != original:
                raise RuntimeError('Completed sample checkpoint changed: ' + sample)
        peak = self.plan.get('resource_probe_peak_rss_bytes', 0)
        self.plan.update(state='active', transition_utc=now(), error=None,
                         estimated_worker_peak_bytes=max(5 * GIB, int(peak * 1.15)),
                         takeover_receipt={'path':str(receipt_path),
                           'sha256':hashlib.sha256(receipt_path.read_bytes()).hexdigest()},
                         original_serial_stopped=True,
                         interrupted_sample='L4', reused_trimmed_sample='L4')
        self.manifest['parallel_execution'] = self.plan.copy()
        self.manifest['state'] = 'running_preprocessing'
        self.manifest.pop('error', None)
        actual_driver = Path(__file__).resolve()
        self.manifest['source_snapshots'][actual_driver.name] = {
            'path':str(actual_driver), 'sha256':hashlib.sha256(actual_driver.read_bytes()).hexdigest()}
        save(self.plan_path, self.plan)
        self.persist()
        self.status('parallel_starting')
        print(now(), 'verified stopped serial queue; two-worker takeover', flush=True)

    def persist(self):
        with self.lock:
            self.manifest['updated_utc'] = now()
            save(self.manifest_path, self.manifest)

    def status(self, state='parallel_libraries', **extra):
        with self.lock:
            save(self.run / 'STATUS.json', {'updated_utc': now(), 'pid': os.getpid(), 'state': state,
                 'execution_mode': 'two_library_workers', 'workers': self.worker_stages.copy(),
                 'completed_samples': list(self.manifest.get('completed_samples', {})), **extra})

    def execute(self):
        try:
            self.handoff()
            spec = importlib.util.spec_from_file_location('worker_driver', self.run / 'worker_driver.py')
            base = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(base)
            coordinator = self

            class Worker(base.Driver):
                def __init__(self, sample):
                    self.root, self.run = coordinator.root, coordinator.run
                    self.refs = self.root / 'references'
                    self.manifest_path = coordinator.manifest_path
                    self.manifest = coordinator.manifest
                    self.mutation_lock = coordinator.lock
                    self.env = os.environ.copy()
                    self.env['PYTHONNOUSERSITE'] = '1'
                    self.env['PATH'] = os.environ.get('PATH', '')
                    self.current_sample = sample
                    self.alignment_threads, self.sort_threads = 3, 1
                    self.trim_cores, self.qc_threads = 1, 1

                def persist(self):
                    coordinator.persist()

                def status(self, stage, **extra):
                    with coordinator.lock:
                        coordinator.worker_stages[self.current_sample] = {'stage': stage, 'updated_utc': now(), **extra}
                        coordinator.status()
                    print(now(), self.current_sample, stage, flush=True)

                def command(self, args, log, stdout_file=None):
                    receipt = json.loads((self.run / 'parallel_takeover_QC.json').read_text())
                    if self.current_sample == 'L4' and str(args[0]) == base.TOOLS['trim']:
                        def without_cores(values):
                            values = list(map(str, values))
                            pos = values.index('--cores')
                            return values[:pos] + values[pos+2:]
                        if without_cores(args) != without_cores(receipt['L4_trim_command']['args']):
                            raise RuntimeError('L4 trim reuse command differs from verified original.')
                        for path in receipt['L4_trim_files']:
                            expected = receipt['verified_files'][path]
                            stat = Path(path).stat()
                            if stat.st_size != expected['bytes'] or stat.st_mtime_ns != expected['mtime_ns']:
                                raise RuntimeError('Verified L4 trimmed file changed: ' + path)
                        with coordinator.lock:
                            self.manifest.setdefault('reused_commands', []).append({
                                'sample':'L4', 'stage':'trimming', 'utc':now(),
                                'source_command':receipt['L4_trim_command'],
                                'validated_receipt':str(self.run / 'parallel_takeover_QC.json')})
                            self.persist()
                        self.status('trimming_reused')
                        return
                    super().command(args, log, stdout_file)

                def fastqc(self, paths, sample, phase):
                    target = self.run / 'qc' / phase
                    expected = [target / (path.name.removesuffix('.gz').removesuffix('.fq') + '_fastqc.zip') for path in paths]
                    valid = False
                    try:
                        valid = bool(expected)
                        for path in expected:
                            with zipfile.ZipFile(path) as archive:
                                valid = valid and archive.testzip() is None
                    except (OSError, zipfile.BadZipFile):
                        pass
                    if valid:
                        self.status('fastqc_' + phase + '_reused')
                        return
                    super().fastqc(paths, sample, phase)

            rows = {row['sample_id']: row for row in self.manifest['sample_metadata']}
            pending = [s for s in self.manifest['queue'] if s not in self.manifest['completed_samples']]
            index = Path(self.manifest['reference']['index_prefix'])
            failures = []
            with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
                active = {}
                while pending or active:
                    for future in list(active):
                        if future.done():
                            sample = active.pop(future)
                            try:
                                future.result()
                            except BaseException as exc:
                                failures.append({'sample': sample, 'error': repr(exc)})
                                print(now(), 'worker failed', sample, repr(exc), flush=True)
                    if failures:
                        pending.clear()  
                    if pending and len(active) < 2:
                        free = available_memory()
                        proc = processes()
                        resident = sum(proc[n]['rss'] for n in descendants(os.getpid(), proc))
                        required = (len(active) + 1) * self.plan['estimated_worker_peak_bytes'] + 4 * GIB
                        if free + resident >= required and shutil.disk_usage(self.run).free >= 60 * GIB:
                            sample = pending.pop(0)
                            active[pool.submit(Worker(sample).sample, rows[sample], index)] = sample
                        else:
                            self.status('waiting_for_worker_resources', free_memory_gib=round(free / GIB, 2),
                                        required_capacity_gib=round(required / GIB, 2))
                    time.sleep(1)
            if failures:
                raise RuntimeError(str(failures))
            worker = Worker('summary')
            self.status('multiqc_summary')
            worker.command([base.TOOLS['multiqc'], str(self.run / 'qc'), str(self.run / 'logs'), str(self.run / 'trimmed'),
                '--outdir', str(self.run / 'qc'), '--filename', 'multiqc_report', '--force'], self.run / 'logs/multiqc.log')
            self.manifest['finished_utc'] = now()
            self.manifest['state'] = 'complete_preprocessing'
            self.persist()
            self.status('complete_preprocessing')
        except BaseException as exc:
            self.plan.update(state='failed', error=repr(exc), updated_utc=now())
            save(self.plan_path, self.plan)
            if self.manifest is not None:
                self.manifest.update(state='failed', error=repr(exc))
                self.persist()
                self.status('failed', error=repr(exc))
            print(now(), 'parallel controller failed', repr(exc), flush=True)
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-root', type=Path, required=True)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    lock = (args.run / '.parallel_controller.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    Coordinator(args.package_root, args.run).execute()
