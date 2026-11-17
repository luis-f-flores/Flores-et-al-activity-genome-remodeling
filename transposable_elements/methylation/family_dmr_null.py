"""Compare TE-family DMR overlap with chromosome-matched shuffles in testable WGBS intervals."""
from pathlib import Path
import argparse, concurrent.futures as cf, csv, json, subprocess, tempfile, time
import numpy as np
import pandas as pd
TIMEPOINTS=('1hr','6hr','24hr')
def run(cmd, *, stdout=None):
    result = subprocess.run(cmd, text=True, stdout=stdout, stderr=subprocess.PIPE)
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}): {' '.join(map(str, cmd))}\n{result.stderr}")
    return result.stderr

def bed_count(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith('#'))

def read_chrom_sizes(path: Path):
    order, sizes = [], {}
    with path.open() as handle:
        for line in handle:
            chrom, size = line.rstrip().split('\t')[:2]
            order.append(chrom)
            sizes[chrom] = int(size)
    return order, sizes

def subtract_intervals(blocks, exclusions):
    """Subtract sorted half-open exclusions from one chromosome's blocks."""
    out = []
    j = 0
    for start, end in blocks:
        while j < len(exclusions) and exclusions[j][1] <= start:
            j += 1
        pieces = [(start, end)]
        k = j
        while k < len(exclusions) and exclusions[k][0] < end:
            es, ee = exclusions[k]
            new = []
            for ps, pe in pieces:
                if ee <= ps or es >= pe:
                    new.append((ps, pe))
                else:
                    if ps < es:
                        new.append((ps, es))
                    if ee < pe:
                        new.append((ee, pe))
            pieces = new
            k += 1
        out.extend((s, e) for s, e in pieces if e > s)
    return out

def load_exclusions(path: Path, chroms):
    values = {chrom: [] for chrom in chroms}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            chrom, start, end = line.rstrip().split('\t')[:3]
            if chrom in values:
                values[chrom].append((int(start), int(end)))
    for chrom in values:
        values[chrom].sort()
    return values

def merge_callable_windows(source: Path, destination: Path, chrom_order, chrom_sizes, exclusions):
    """Stream the assessed 1-based closed windows and write masked BED blocks."""
    merged = {chrom: [] for chrom in chrom_order}
    current = {}
    n_windows = 0
    with source.open() as handle:
        header = handle.readline().rstrip().split('\t')
        idx = {name: header.index(name) for name in ('chr', 'start', 'end')}
        for line in handle:
            fields = line.rstrip().split('\t')
            chrom = fields[idx['chr']]
            if chrom not in chrom_sizes:
                continue
            start = int(fields[idx['start']]) - 1
            end = int(fields[idx['end']])
            n_windows += 1
            if chrom not in current:
                current[chrom] = [start, end]
            elif start <= current[chrom][1]:
                if end > current[chrom][1]:
                    current[chrom][1] = end
            else:
                merged[chrom].append(tuple(current[chrom]))
                current[chrom] = [start, end]
    for chrom, interval in current.items():
        merged[chrom].append(tuple(interval))

    masked = {chrom: subtract_intervals(merged[chrom], exclusions[chrom]) for chrom in chrom_order}
    with destination.open('w') as out:
        for chrom in chrom_order:
            for start, end in masked[chrom]:
                out.write(f'{chrom}\t{start}\t{end}\n')
    return {
        'assessed_windows': n_windows,
        'merged_blocks_before_mask': sum(len(v) for v in merged.values()),
        'blocks_after_mask': sum(len(v) for v in masked.values()),
        'bp_after_mask': sum(e - s for v in masked.values() for s, e in v),
    }

def build_dmr_bed(sources: dict, tp: str, raw_bed: Path, eligible_bed: Path,
                  chrom_sizes_path: Path, exclusions_path: Path):
    _, chrom_sizes = read_chrom_sizes(chrom_sizes_path)
    inputs = sources[tp]
    n = 0
    with raw_bed.open('w') as out:
        for direction in ('hyper', 'hypo'):
            frame = pd.read_csv(inputs[direction], sep='\t', usecols=['chr', 'start', 'end'])
            frame = frame[frame['chr'].isin(chrom_sizes)]
            for row in frame.itertuples(index=False):
                start0 = int(row.start) - 1
                end0 = int(row.end)
                if not (0 <= start0 < end0 <= chrom_sizes[row.chr]):
                    raise ValueError(f'Invalid DMR coordinates: {row}')
                out.write(f'{row.chr}\t{start0}\t{end0}\t{tp}_{direction}_{n}\t0\t{direction}\n')
                n += 1
    with eligible_bed.open('w') as out:
        run(['bedtools', 'intersect', '-v', '-a', str(raw_bed), '-b', str(exclusions_path)], stdout=out)
    return {'raw_canonical_dmrs': n, 'eligible_masked_dmrs': bed_count(eligible_bed)}

def top_families(table: Path):
    frame = pd.read_csv(table)
    scores = (frame.groupby('TE_family', as_index=False)['pct_copies_overlapping']
              .sum().rename(columns={'pct_copies_overlapping': 'panel_B_total_pct'}))
    scores = scores.sort_values(['panel_B_total_pct', 'TE_family'], ascending=[False, True])
    cutoff = scores.iloc[min(19, len(scores)-1)]['panel_B_total_pct']
    chosen = scores[scores['panel_B_total_pct'] >= cutoff].copy()
    return chosen

def build_top_te_bed(te_meta: Path, selected: pd.DataFrame, destination: Path):
    wanted = set(selected['TE_family'])
    totals = {}
    kept = 0
    with te_meta.open() as handle, destination.open('w') as out:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            fam = row['family']
            totals[fam] = totals.get(fam, 0) + 1
            if fam in wanted:
                out.write(f"{row['chr']}\t{row['start']}\t{row['end']}\t{row['copyID']}\t{fam}\n")
                kept += 1
    return totals, kept

def observed_counts(top_te_bed: Path, dmr_beds, output_bed: Path):
    with output_bed.open('w') as out:
        run(['bedtools', 'intersect', '-u', '-a', str(top_te_bed), '-b', *map(str, dmr_beds)], stdout=out)
    counts = {}
    with output_bed.open() as handle:
        for line in handle:
            fam = line.rstrip().split('\t')[4]
            counts[fam] = counts.get(fam, 0) + 1
    return counts

def one_iteration(task):
    model, iteration, seed, work, chrom_sizes, top_te, dmr_beds, shuffle_exclusions = task
    started = time.time()
    with tempfile.TemporaryDirectory(prefix=f'te_dmr_{model}_{iteration}_') as td:
        td = Path(td)
        shuffled = []
        input_counts = []
        output_counts = []
        for j, tp in enumerate(TIMEPOINTS):
            out = td / f'{tp}.bed'
            cmd = ['bedtools', 'shuffle', '-i', str(dmr_beds[tp]), '-g', str(chrom_sizes),
                   '-excl', str(shuffle_exclusions[tp]),
                   '-chrom', '-noOverlapping', '-maxTries', '1000000',
                   '-seed', str(seed + iteration * 101 + j)]
            with out.open('w') as handle:
                run(cmd, stdout=handle)
            input_counts.append(bed_count(dmr_beds[tp]))
            output_counts.append(bed_count(out))
            if input_counts[-1] != output_counts[-1]:
                raise RuntimeError(f'{model} iteration {iteration} {tp}: lost shuffled intervals')
            shuffled.append(out)
        hits = td / 'hits.bed'
        with hits.open('w') as handle:
            run(['bedtools', 'intersect', '-u', '-a', str(top_te), '-b', *map(str, shuffled)], stdout=handle)
        counts = {}
        with hits.open() as handle:
            for line in handle:
                fam = line.rstrip().split('\t')[4]
                counts[fam] = counts.get(fam, 0) + 1
    return model, iteration, counts, time.time() - started

def percentile(values, q):
    return float(np.quantile(np.asarray(values, dtype=float), q))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--settings',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--prepare-only',action='store_true')
    a=p.parse_args();settings=json.loads(a.settings.read_text())
    def source(v):
        path=Path(v)
        return path if path.is_absolute() else (a.settings.parent/path).resolve()
    sources={tp:{k:source(v) for k,v in settings['dmrs'][tp].items()} for tp in TIMEPOINTS}
    chrom_sizes=source(settings['chrom_sizes']);exclusions_path=source(settings['exclusions']);te_meta=source(settings['te_meta'])
    family_counts=source(settings['family_overlap']);out=a.output
    work=out/'work';tables=out/'tables';work.mkdir(parents=True,exist_ok=True);tables.mkdir(parents=True,exist_ok=True)
    chrom_order,chrom_size_dict=read_chrom_sizes(chrom_sizes);exclusions=load_exclusions(exclusions_path,chrom_order)
    selected=top_families(family_counts);selected.to_csv(tables/'top_family_selection.tsv',sep='\t',index=False)
    top_te=work/'top_TE_copies.bed';totals,top_copy_count=build_top_te_bed(te_meta,selected,top_te)
    dmr_beds={};shuffle_exclusions={};support=[]
    for tp in TIMEPOINTS:
        raw=work/f'dmr_{tp}_canonical.bed';eligible=work/f'dmr_{tp}_eligible.bed'
        dmr_stats=build_dmr_bed(sources,tp,raw,eligible,chrom_sizes,exclusions_path);dmr_beds[tp]=eligible
        callable_bed=work/f'callable_{tp}_masked.bed'
        call_stats=merge_callable_windows(sources[tp]['assessed'],callable_bed,chrom_order,chrom_size_dict,exclusions)
        excluded=work/f'callable_{tp}_excluded.bed'
        with excluded.open('w') as handle:run(['bedtools','complement','-i',str(callable_bed),'-g',str(chrom_sizes)],stdout=handle)
        shuffle_exclusions[tp]=excluded
        contained=work/f'dmr_{tp}_fully_callable.bed'
        with contained.open('w') as handle:run(['bedtools','intersect','-u','-f','1.0','-a',str(eligible),'-b',str(callable_bed)],stdout=handle)
        if bed_count(contained)!=bed_count(eligible):raise RuntimeError(f'{tp}: DMRs must be fully contained in the assessed space.')
        support.append(dict(timepoint=tp,**dmr_stats,**call_stats))
    pd.DataFrame(support).to_csv(tables/'input_support.tsv',sep='\t',index=False)
    observed=observed_counts(top_te,[dmr_beds[t] for t in TIMEPOINTS],work/'observed_TE_hits.bed')
    obs=pd.DataFrame([dict(TE_family=r.TE_family,total_catalog_copies=totals[r.TE_family],observed_hit_copies=observed.get(r.TE_family,0),observed_pct=100*observed.get(r.TE_family,0)/totals[r.TE_family]) for r in selected.itertuples(index=False)])
    obs.to_csv(tables/'observed_family_overlap.tsv',sep='\t',index=False)
    if a.prepare_only:return
    iterations=1000;seed=20260925
    tasks=[('callable',i,seed,work,chrom_sizes,top_te,dmr_beds,shuffle_exclusions) for i in range(1,iterations+1)]
    values=[]
    with cf.ProcessPoolExecutor(max_workers=a.workers) as pool:
        for done,(model,iteration,counts,seconds) in enumerate(pool.map(one_iteration,tasks),1):
            for fam in selected.TE_family:
                count=counts.get(fam,0);values.append(dict(null_model=model,iteration=iteration,TE_family=fam,hit_copies=count,pct_copies=100*count/totals[fam]))
            if done%10==0 or done==iterations:print(f'Shuffles completed: {done}/{iterations}',flush=True)
    null=pd.DataFrame(values);null.to_csv(tables/'family_null_iterations.tsv.gz',sep='\t',index=False)
    rows=[]
    for r in obs.itertuples(index=False):
        v=null.loc[null.TE_family==r.TE_family,'pct_copies'].to_numpy();sd=float(np.std(v,ddof=1))
        rows.append(dict(null_model='callable',TE_family=r.TE_family,total_catalog_copies=r.total_catalog_copies,observed_hit_copies=r.observed_hit_copies,observed_pct=r.observed_pct,
          null_mean_pct=float(np.mean(v)),null_sd_pct=sd,null_q025_pct=percentile(v,.025),null_q975_pct=percentile(v,.975),
          observed_minus_null_pct_points=r.observed_pct-float(np.mean(v)),null_z=(r.observed_pct-float(np.mean(v)))/sd if sd>0 else float('nan'),
          empirical_p_enrichment_selection_conditioned=(1+int(np.sum(v>=r.observed_pct)))/(len(v)+1),
          empirical_p_depletion_selection_conditioned=(1+int(np.sum(v<=r.observed_pct)))/(len(v)+1)))
    pd.DataFrame(rows).to_csv(tables/'34_te_family_null_summary.tsv',sep='\t',index=False)

if __name__=='__main__':main()
