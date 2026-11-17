
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
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def record(path, full=True):
    p = Path(path).resolve()
    s = p.stat()
    return dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns,
                sha256=sha256(p) if full else None)


def save(path, data):
    p = Path(path)
    tmp = p.with_suffix(p.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    tmp.replace(p)


def table(path, header, rows):
    with Path(path).open('w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(header)
        w.writerows(rows)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def alias(chrom):
    """Only known mm10 primary/accession naming differences; not liftover."""
    if chrom in ('M', 'MT', 'chrM'):
        return 'MT'
    m = re.fullmatch(r'(?:chr)?([1-9]|1[0-9]|X|Y)', chrom)
    if m:
        return m.group(1)
    m = re.search(r'(GL\d+|JH\d+)(?:\.\d+)?', chrom)
    return m.group(1) if m else chrom


def closed_to_bed(start, end):
    require(1 <= start <= end, 'Invalid 1-based closed interval')
    return start - 1, end


def bed_to_closed(start, end):
    require(0 <= start < end, 'Invalid 0-based half-open interval')
    return start + 1, end


def sample_rows(raw_dir, bam_dir):
    rows = []
    for tp in (0, 1, 6, 24):
        for treatment in (['baseline'] if tp == 0 else ['M', 'K']):
            for rep in 'ABC':
                name = f'{rep}{tp}' + ('' if tp == 0 else treatment)
                rows.append(dict(sample=name, replicate=rep, legacy_timepoint=str(tp),
                                 timepoint='recovery' if tp == 24 else f'{tp}hr',
                                 condition={'M': 'mannitol', 'K': 'KCl'}.get(treatment, treatment),
                                 washout_to_harvest='12-16hr' if tp == 24 else 'NA',
                                 fastq1=str(Path(raw_dir) / (name + '_1.fq.gz')),
                                 fastq2=str(Path(raw_dir) / (name + '_2.fq.gz')),
                                 original_bam=str(Path(bam_dir) / (name + '.bam'))))
    return rows


def version(program):
    r = subprocess.run([program, '--version'], capture_output=True, text=True,
                       errors='replace', check=True)
    return (r.stdout + r.stderr).strip()


def verify_records(records):
    for r in records:
        p = Path(r['path'])
        require(p.is_file() and p.stat().st_size == r['bytes'], f'Changed/missing input: {p}')
        if r.get('sha256'):
            require(sha256(p) == r['sha256'], f'Hash mismatch: {p}')
        else:
            require(p.stat().st_mtime_ns == r['mtime_ns'], f'Changed input timestamp: {p}')


def disk_check(path, reserve_gib):
    free = shutil.disk_usage(path).free / 2**30
    require(free >= reserve_gib, f'Free space {free:.1f} GiB < reserve {reserve_gib} GiB')
    return free
