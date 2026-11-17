'Prepare inputs for the ChIP mapping driver.'
from pathlib import Path
import argparse,csv,json
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--raw',type=Path,required=True,help='Directory containing the delivered FASTQs, possibly in sample subdirectories.')
p.add_argument('--mouse-fasta',type=Path,required=True)
p.add_argument('--run',type=Path,required=True)
a=p.parse_args();a.run=a.run.resolve()
if a.run.exists():raise FileExistsError('Choose a new run directory.')
with Path(__file__).with_name('raw_reads.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
files=[];samples=[]
for sid in dict.fromkeys(r['sample_id'] for r in rows):
    pair=sorted([r for r in rows if r['sample_id']==sid],key=lambda r:r['filename'])
    assert len(pair)==2
    paths=[]
    for r in pair:
        found=list(a.raw.rglob(r['filename']))
        if len(found)!=1:raise FileNotFoundError(r['filename']+' must occur exactly once.')
        q=found[0].resolve();paths.append(str(q))
        files.append(dict(sample_id=sid,local_path=str(q),expected_bytes=int(r['expected_bytes']),expected_md5=r['expected_md5']))
    r=pair[0];samples.append(dict(sample_id=sid,condition=r['condition'],bio_rep=r['biological_replicate'],assay=r['assay'],r1_path=paths[0],r2_path=paths[1]))
assert a.mouse_fasta.is_file() and Path(str(a.mouse_fasta)+'.fai').is_file()
for name in ['bam','trimmed','qc','logs']:(a.run/name).mkdir(parents=True)
d=dict(settings=dict(mouse_fasta=str(a.mouse_fasta.resolve()),mouse_fai=str(a.mouse_fasta.resolve())+'.fai'),sample_metadata=samples,raw_files=files,queue=[s['sample_id'] for s in samples],completed_samples={},commands=[])
(a.run/'RUN_MANIFEST.json').write_text(json.dumps(d,indent=2)+'\n')
print('Prepared '+str(a.run)+'; no processing launched.')
