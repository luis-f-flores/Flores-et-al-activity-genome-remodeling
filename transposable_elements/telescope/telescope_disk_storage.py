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

import ast
import gc
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import pickle
import signal
import sys
import textwrap
import time

from telescope_common import now, record, require, save, verify_records

PACKAGE = GB_WORK_ROOT / 'rna_package'
MODEL_SHA256 = '13d7031334ff777f9a0fed592c0b82c1dda1da4466924a31ae71cec974f6ec31'


def module(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent/f'{name}.py')
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


class DiskMappings:
    instances = []

    def __init__(self, opts):
        self.path = Path(opts.outfile_path('mappings.tmp'))
        self.stream = self.path.open('xb', buffering=1024*1024)
        self.count = self.bytes = self.read_count = 0
        self.iterated = False
        self.instances.append(self)

    def append(self, item):
        require(not self.iterated, 'Cannot append after reading mapping spool')
        
        data = pickle.dumps(item, protocol=4)
        self.stream.write(data)
        self.bytes += len(data)
        self.count += 1

    def __iter__(self):
        require(not self.iterated, 'Mapping spool must be consumed exactly once')
        self.iterated = True
        self.stream.flush()
        with self.path.open('rb', buffering=1024*1024) as source:
            for _ in range(self.count):
                value = pickle.load(source)
                self.read_count += 1
                yield value
            require(source.read(1) == b'', 'Extra bytes in mapping spool')
        require(self.read_count == self.count, 'Incomplete mapping spool read')

    def close(self):
        self.stream.close()
        
        
        self.path.unlink(missing_ok=True)


def install_disk_storage():
    import telescope.utils.model as model
    require(record(model.__file__)['sha256'] == MODEL_SHA256, 'Telescope model source checksum does not match the required version')
    original = model.Telescope._load_sequential
    tree = ast.parse(textwrap.dedent(inspect.getsource(original)))
    replacements = 0
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and len(node.targets)==1 and
                isinstance(node.targets[0], ast.Name) and node.targets[0].id=='_mappings' and
                isinstance(node.value, ast.List) and not node.value.elts):
            node.value = ast.Call(func=ast.Name(id='_DiskBackedMappings',ctx=ast.Load()),
                                 args=[ast.Attribute(value=ast.Name(id='self',ctx=ast.Load()),
                                                     attr='opts',ctx=ast.Load())],keywords=[])
            replacements += 1
    require(replacements == 1, 'Expected exactly one mapping-list allocation')
    ast.fix_missing_locations(tree)
    
    
    model.__dict__['_DiskBackedMappings'] = DiskMappings
    local = {}
    exec(compile(tree, str(Path(__file__))+'::<storage-only adaptation>', 'exec'),model.__dict__,local)
    model.Telescope._load_sequential = local['_load_sequential']
    return model, original


def assign():
    model, original = install_disk_storage()
    cached = module('telescope_annotation_cache')
    outdir = Path(sys.argv[sys.argv.index('--outdir')+1])
    tag = sys.argv[sys.argv.index('--exp_tag')+1]
    output = outdir/(tag+'.loading_storage.json')
    require(not output.exists(), 'Loading-storage QC already exists; use guarded recovery')
    def interrupted(signum, frame):
        raise InterruptedError('Loading interrupted by signal '+str(signum))
    for signum in (signal.SIGTERM,signal.SIGINT,signal.SIGHUP):
        signal.signal(signum,interrupted)
    started=time.monotonic();status='failed'
    try:
        cached.run()
        require(all(x.count == x.read_count for x in DiskMappings.instances), 'Incomplete mapping consumption')
        status='complete'
    finally:
        stats=[dict(mapping_tuples=x.count,consumed_tuples=x.read_count,
                    temporary_bytes=x.bytes,path=str(x.path)) for x in DiskMappings.instances]
        for x in DiskMappings.instances:x.close()
        save(output,dict(status=status,finished_utc=now(),source=record(__file__),
                         elapsed_seconds=time.monotonic()-started,storage=stats,
                         temporary_files_removed=True,
                         adaptation='Stream mapping-list entries from disk in input order for Telescope matrix construction and assignment.'))
        model.Telescope._load_sequential=original



if __name__ == "__main__":
    assign()
