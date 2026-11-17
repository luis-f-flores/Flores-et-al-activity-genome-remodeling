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

import collections
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone
import pysam

ROOT = Path(gb_path('chipseq'))
RUN = ROOT / 'runs/2026-10-05_trim_map'
THRESHOLDS = (0, 1, 10, 20, 30)


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def median(histogram):
    total = sum(histogram.values())
    cumulative = 0
    for value, n in sorted(histogram.items()):
        cumulative += n
        if cumulative >= (total + 1) // 2:
            return value
    return None


def inspect(sample):
    sid = sample['sample_id']
    bam_path = Path(sample['bam_path'])
    expected_rows = list(csv.DictReader((RUN / 'qc/library_qc.tsv').open(), delimiter='\t'))
    expected = next(row for row in expected_rows if row['sample_id'] == sid)
    counters = collections.Counter()
    mapq = collections.Counter()
    pair_mapq = collections.Counter()
    insert_sizes = collections.Counter()
    pair_counts = collections.Counter()
    contig_counts = collections.Counter()
    before = (bam_path.stat().st_size, bam_path.stat().st_mtime_ns)
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        fly = [name for name in bam.references if name.startswith('dm6_')]
        fly_ids = {bam.get_tid(name) for name in fly}
        for contig in fly:
            for read in bam.fetch(contig):
                if read.is_unmapped:
                    continue
                counters['fly_mapped_read_ends'] += 1
                contig_counts[contig] += 1
                if read.is_duplicate:
                    counters['duplicate_read_ends'] += 1
                if read.is_secondary or read.is_supplementary:
                    counters['nonprimary_read_ends'] += 1
                    continue
                counters['primary_read_ends'] += 1
                mapq[read.mapping_quality] += 1
                if not read.is_duplicate:
                    counters['nonduplicate_primary_read_ends'] += 1
                if read.is_qcfail:
                    counters['qc_failed_primary_read_ends'] += 1
                if read.mate_is_unmapped:
                    counters['mate_unmapped_primary_read_ends'] += 1
                elif read.next_reference_id not in fly_ids:
                    counters['mate_mouse_primary_read_ends'] += 1
                elif read.is_proper_pair:
                    counters['proper_fly_pair_primary_read_ends'] += 1
                else:
                    counters['nonproper_fly_pair_primary_read_ends'] += 1
                proper = (read.is_paired and read.is_proper_pair and not read.mate_is_unmapped
                          and read.next_reference_id == read.reference_id and not read.is_qcfail)
                if not proper:
                    continue
                if not read.has_tag('MQ'):
                    counters['proper_pair_missing_mate_mapq_read_ends'] += 1
                    continue
                if not read.template_length:
                    counters['proper_pair_zero_tlen_read_ends'] += 1
                    continue
                minimum_q = min(read.mapping_quality, read.get_tag('MQ'))
                mate = 1 if read.is_read1 else 2 if read.is_read2 else 0
                assert mate, f'{sid}: primary proper pair has no mate bit'
                if mate == 1:
                    pair_mapq[minimum_q] += 1
                    insert_sizes[abs(read.template_length)] += 1
                groups = ['all_fly']
                if contig != 'dm6_chrM':
                    groups.append('nonmitochondrial_fly')
                modes = ['duplicates_retained']
                if not read.is_duplicate:
                    modes.append('duplicates_excluded')
                for group in groups:
                    for mode in modes:
                        for cutoff in THRESHOLDS:
                            if minimum_q >= cutoff:
                                pair_counts[group, mode, cutoff, mate] += 1
    assert before == (bam_path.stat().st_size, bam_path.stat().st_mtime_ns)
    total = counters['fly_mapped_read_ends']
    assert total == int(expected['fly_mapped_read_ends']), (sid, total, expected)
    assert counters['proper_pair_missing_mate_mapq_read_ends'] == 0
    support = []
    for group in ['all_fly', 'nonmitochondrial_fly']:
        for mode in ['duplicates_retained', 'duplicates_excluded']:
            for cutoff in THRESHOLDS:
                n1 = pair_counts[group, mode, cutoff, 1]
                n2 = pair_counts[group, mode, cutoff, 2]
                assert n1 == n2, (sid, group, mode, cutoff, n1, n2)
                support.append({'sample_id': sid, 'condition': sample['condition'],
                                'bio_rep': sample['bio_rep'], 'contigs': group,
                                'duplicate_handling': mode, 'minimum_MAPQ_both_mates': cutoff,
                                'proper_paired_fragments': n1})
    
    comparisons = []
    for args, key in [(['/usr/bin/samtools', 'view', '-c', '-F', '4', str(bam_path)] + fly,
                       'fly_mapped_read_ends'),
                      (['/usr/bin/samtools', 'view', '-c', '-F', '4', '-f', '1024', str(bam_path)] + fly,
                       'duplicate_read_ends')]:
        count = int(subprocess.run(args, capture_output=True, text=True, check=True).stdout)
        assert count == counters[key], (sid, key, count, counters[key])
        comparisons.append({'args': args, 'counter': key, 'count': count, 'matched': True})
    summary = {'sample_id': sid, 'condition': sample['condition'], 'bio_rep': sample['bio_rep'],
               'fly_mapped_read_ends': total,
               'fly_duplicate_read_ends_pct': round(100 * counters['duplicate_read_ends'] / total, 4),
               'fly_primary_MAPQ0_pct': round(100 * mapq[0] / counters['primary_read_ends'], 4),
               'fly_primary_MAPQ_ge10_pct': round(100 * sum(n for q, n in mapq.items() if q >= 10) / counters['primary_read_ends'], 4),
               'fly_primary_MAPQ_ge30_pct': round(100 * sum(n for q, n in mapq.items() if q >= 30) / counters['primary_read_ends'], 4),
               'fly_proper_pair_primary_read_ends_pct': round(100 * counters['proper_fly_pair_primary_read_ends'] / counters['primary_read_ends'], 4),
               'fly_mate_mouse_primary_read_ends': counters['mate_mouse_primary_read_ends'],
               'fly_mitochondrial_read_ends_pct': round(100 * contig_counts['dm6_chrM'] / total, 4),
               'fly_primary_MAPQ_median': median(mapq), 'fly_pair_min_MAPQ_median': median(pair_mapq),
               'fly_proper_pair_insert_median_bp': median(insert_sizes),
               'fly_nonduplicate_proper_pairs_MAPQ0': pair_counts['all_fly', 'duplicates_excluded', 0, 1],
               'fly_nonduplicate_proper_pairs_MAPQ10': pair_counts['all_fly', 'duplicates_excluded', 10, 1],
               'fly_nonduplicate_proper_pairs_MAPQ20': pair_counts['all_fly', 'duplicates_excluded', 20, 1],
               'fly_nonduplicate_proper_pairs_MAPQ30': pair_counts['all_fly', 'duplicates_excluded', 30, 1],
               'fly_nonmitochondrial_nonduplicate_proper_pairs_MAPQ10': pair_counts['nonmitochondrial_fly', 'duplicates_excluded', 10, 1]}
    detail = {'sample_id': sid, 'bam_path': str(bam_path), 'bam_size': before[0],
              'bam_unchanged_size_mtime': True, 'counters': dict(counters),
              'read_end_MAPQ_histogram': dict(mapq), 'proper_pair_min_MAPQ_histogram': dict(pair_mapq),
              'proper_pair_insert_histogram_bp': dict(insert_sizes), 'mapped_read_ends_by_contig': dict(contig_counts),
              'independent_samtools_count_checks': comparisons,
              'read1_read2_pair_counts_match_all_thresholds': True}
    return summary, support, detail


def main():
    result_path = RUN / 'qc/15_FLY_QC_RESULT.json'
    assert not result_path.exists(), 'Inspect existing QC before repeating this diagnostic.'
    samples = [s for s in csv.DictReader((ROOT / 'samples.tsv').open(), delimiter='\t') if s['assay'] == 'IP']
    assert {s['sample_id'] for s in samples} == {f'L{i}' for i in range(9,17)}
    source = Path(__file__).resolve()
    record = {'started_utc': now(), 'state': 'running', 'source': str(source),
              'source_sha256': sha(source), 
              'versions': {'pysam': pysam.__version__, 'pysam_samtools': pysam.__samtools_version__},
              'workers': 2, 'MAPQ_thresholds_for_diagnostics_only': list(THRESHOLDS),
              'paired_fragment_definition': 'One primary read1 of a QC-passing proper pair on the same fly contig; both mates meet MAPQ cutoff via fixmate MQ tag; nonzero TLEN. Duplicate-retained/excluded and all/nonmitochondrial fly are reported separately.',
              
              
              'input_evidence': []}
    checkpoints = {}
    for m in [RUN/'RUN_MANIFEST.json', RUN/'beast_transfer/RUN_MANIFEST.json']:
        checkpoints.update(json.loads(m.read_text())['completed_samples'])
    for sample in samples:
        bam_path=Path(sample['bam_path']);cp=checkpoints[sample['sample_id']]
        cp_path=next(path for path in cp['outputs'] if path.endswith('.marked.bam'))
        evidence=cp['outputs'][cp_path]
        assert bam_path.stat().st_size == evidence['bytes']
        record['input_evidence'].append({'sample_id':sample['sample_id'],'local_bam':str(bam_path),
                    'previously_verified_bam_sha256':evidence['sha256'],'previously_verified_bam_bytes':evidence['bytes'],
                    'idxstats_sha256':sha(RUN/'qc'/f"{sample['sample_id']}.idxstats.txt")})
    record['sample_key_sha256'] = sha(ROOT/'samples.tsv')
    result_path.write_text(json.dumps(record,indent=2)+'\n')
    summaries,support,details=[],[],[]
    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=2) as pool:
            jobs={pool.submit(inspect,s):s['sample_id'] for s in samples}
            for job in concurrent.futures.as_completed(jobs):
                summary,pairs,detail=job.result()
                summaries.append(summary);support.extend(pairs);details.append(detail)
                print(json.dumps(summary),flush=True)
        summaries.sort(key=lambda x:int(x['sample_id'][1:]));support.sort(key=lambda x:(int(x['sample_id'][1:]),x['contigs'],x['duplicate_handling'],x['minimum_MAPQ_both_mates']))
        details.sort(key=lambda x:int(x['sample_id'][1:]))
        summary_path=RUN/'qc/fly_qc.tsv';pairs_path=RUN/'qc/fly_pair_support.tsv'
        write_tsv(summary_path,summaries);write_tsv(pairs_path,support)
        record.update({'state':'complete','completed_utc':now(),'samples':details,
                       'independent_count_checks_passed':16,'all_pair_thresholds_read1_read2_match':True,
                       'outputs':[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)} for p in [summary_path,pairs_path]]})
        result_path.write_text(json.dumps(record,indent=2)+'\n')
        print('COMPLETE: eight-IP read-only fly diagnostics; 16 SAMtools recounts and all mate-count checks pass.',flush=True)
    except Exception as error:
        record.update({'state':'failed','error':str(error),'failed_utc':now()})
        result_path.write_text(json.dumps(record,indent=2)+'\n')
        raise


if __name__=='__main__':
    main()
