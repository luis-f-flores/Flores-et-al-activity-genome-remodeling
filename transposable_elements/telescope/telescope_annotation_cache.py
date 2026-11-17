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
from collections import Counter
import hashlib
import io
import json
import os
from pathlib import Path
import resource
import sys
import time
from telescope_common import record, require, save, now, verify_records
from telescope.utils._annotation_intervaltree import _AnnotationIntervalTree as Annotation
from intervaltree import Interval


def signature(a):
    h=hashlib.sha256()
    for name,features in a.loci.items():
        h.update((name+'\n').encode())
        for f in features:h.update(('\t'.join(f)+'\n').encode())
    for chrom in sorted(a.itree):
        for iv in sorted(a.itree[chrom],key=lambda x:(x.begin,x.end,x.data[a.key])):
            h.update(f'{chrom}\t{iv.begin}\t{iv.end}\t{iv.data[a.key]}\n'.encode())
    return h.hexdigest()


def strict_overlap(a,ref,blocks):
    
    
    result=Counter()
    for start,end in blocks:
        require(end>=start,'Reversed read block')
        if end==start:continue
        query=Interval(start+1,end+1)
        for iv in a.itree[ref].overlap(query):
            result[iv.data[a.key]]+=max(0,min(iv.end,query.end)-max(iv.begin,query.begin))
    return result


def build(gtf,cache):
    require(not cache.exists() and not cache.with_suffix('.json').exists(),'Cache already exists')
    cache.parent.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    with gtf.open() as handle:
        a=Annotation(handle,'locus')
    before=signature(a);loci=len(a.loci)
    a.save(str(cache));del a
    b=Annotation.load(str(cache));require(signature(b)==before,'Serialized annotation differs')
    toy=Annotation(io.StringIO('chrT\tX\texon\t10\t20\t.\t+\t.\tlocus "test";\n'),'locus')
    expected={'full':((9,20),11),'inside':((12,18),6),'right':((20,21),0),'left':((8,9),0)}
    for label,(block,length) in expected.items():
        require(strict_overlap(toy,'chrT',[block]).get('test',0)==length,f'Boundary fixture failed: {label}')
    save(cache.with_suffix('.json'),dict(created_utc=now(),source=record(gtf),cache=record(cache),
        code=record(Path(__file__)),loci=loci,annotation_signature=before,
        serialization_roundtrip_identical=True,boundary_fixture_checks=4,
        elapsed_seconds=time.monotonic()-started,
        peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
    print(json.dumps(dict(status='cache_verified',loci=loci,cache=str(cache)),indent=2),flush=True)


def run():
    import telescope.__main__ as cli
    import telescope.telescope_assign as assign
    if 'assign' not in sys.argv or '--help' in sys.argv:
        return cli.main()
    cache=Path(os.environ['TELESCOPE_ANNOTATION_CACHE'])
    metadata=json.loads(cache.with_suffix('.json').read_text())
    verify_records([metadata['cache'],metadata['source']])
    class CachedAnnotation:
        def __init__(self,gtf,attribute):
            require(attribute=='locus','Unexpected locus key')
            require(Path(gtf).resolve()==Path(metadata['source']['path']).resolve(),'Cache is for another GTF')
            self.annotation=Annotation.load(str(cache));self.loci=self.annotation.loci
        def feature_length(self): return self.annotation.feature_length()
        def intersect_blocks(self,reference,blocks): return self.annotation.intersect_blocks(reference,blocks)
    assign.get_annotation_class=lambda name:CachedAnnotation
    cli.main()


if __name__=='__main__':
    if '--build-cache' in sys.argv:
        p=argparse.ArgumentParser();p.add_argument('--build-cache',type=Path,required=True);p.add_argument('--gtf',type=Path,required=True)
        a=p.parse_args();build(a.gtf,a.build_cache)
    else:run()
