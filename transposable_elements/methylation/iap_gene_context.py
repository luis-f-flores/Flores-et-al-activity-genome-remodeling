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

import bisect, collections, csv, gzip, json, re, resource
from pathlib import Path
from telescope_common import now, record, require, verify_records, table, alias

PACKAGE = GB_WORK_ROOT / 'rna_package'
CHECK=PACKAGE/'runs/2026-09-14_telescope_full6h_benchmark/checks'
REVIEW=PACKAGE/'runs/2026-09-13_iap_locus_validation'
OUT='C6_IAP_gene_context'
def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def overlap(a,b,c,d):return max(0,min(b,d)-max(a,c))
def merged_length(intervals):
    total=0;left=right=None
    for a,b in sorted(intervals):
        if right is None:left,right=a,b
        elif a<=right:right=max(right,b)
        else:total+=right-left;left,right=a,b
    return total+(right-left if right is not None else 0)
def parse(path,loci,lengths):
    targets=collections.defaultdict(list)
    for lid,r in loci.items():targets[r['chr']].append((r['start0'],r['end'],lid))
    for c in targets:targets[c].sort()
    starts={c:[r[0] for r in rr] for c,rr in targets.items()}
    genes={}; exons=collections.defaultdict(lambda:collections.defaultdict(list))
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as f:
        for line in f:
            if line.startswith('#'):continue
            v=line.rstrip().split('\t')
            require(len(v)==9,'Malformed GTF row')
            c=alias(v[0]);feature=v[2]
            if c not in targets or feature not in ('gene','transcript','exon'):continue
            a,b=int(v[3])-1,int(v[4]);strand=v[6]
            require(0<=a<b<=lengths[c] and strand in ('+','-'),'Invalid GTF coordinates')
            hits=[(s,e,lid) for s,e,lid in targets[c][:bisect.bisect_left(starts[c],b)] if e>a] if feature=='exon' else []
            if feature=='exon' and not hits:continue
            attrs=dict(re.findall(r'(\w+) "([^"]*)"',v[8]))
            gid=attrs['gene_id'];key=(c,strand,gid)
            if feature in ('gene','transcript'):
                g=genes.setdefault(key,dict(chr=c,strand=strand,gid=gid,name=attrs.get('gene_name',gid),start0=a,end=b))
                g['start0']=min(g['start0'],a);g['end']=max(g['end'],b)
            else:
                for s,e,lid in hits:exons[lid][gid].append((max(s,a),min(e,b)))
    bychr=collections.defaultdict(list)
    for g in genes.values():bychr[g['chr']].append(g)
    result={}
    for lid,r in loci.items():
        body={g['gid']:dict(g,body=overlap(r['start0'],r['end'],g['start0'],g['end']),exon=merged_length(exons[lid].get(g['gid'],[])))
              for g in bychr[r['chr']] if overlap(r['start0'],r['end'],g['start0'],g['end'])}
        exon_bp=merged_length([i for intervals in exons[lid].values() for i in intervals])
        if exon_bp:category='exon_overlapping'
        elif not body:category='outside_gene_bodies'
        elif any(g['start0']<=r['start0'] and g['end']>=r['end'] for g in body.values()):category='fully_intronic'
        else:category='gene_boundary_overlap'
        result[lid]=dict(category=category,body=body,exon_bp=exon_bp,
            body_genes=','.join(sorted(set(g['name'] for g in body.values()))),
            exon_genes=','.join(sorted(set(g['name'] for g in body.values() if g['exon']))))
    return result
def summarize(rr):
    categories=dict(collections.Counter(r['context_both_annotations'] for r in rr))
    return dict(n=len(rr),categories=categories,loci=[r['locus'] for r in rr])
def main():
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    require(not any((CHECK/(OUT+ext)).exists() for ext in ('.tsv','.json')),'Preserve completed outputs')
    pq=CHECK/'C6_pair_review.json';q=json.loads(pq.read_text());verify_records(q['outputs']+[q['source']])
    gq=REVIEW/'GENE_CONTEXT_QC.json';gqc=json.loads(gq.read_text());verify_records(gqc['inputs']+gqc['outputs']+[x['input'] for x in gqc['annotations'].values()])
    all_rows=rows(CHECK/'C6_pair_IAP_all.tsv')
    pool=[r for r in all_rows if float(r['C6K_final_conf'])>=10 and float(r['C6K_init_best'])>=10]
    scan=json.loads((CHECK/'C6K.support_scan.json').read_text())
    expected=next(t['joint_support_ids'] for t in scan['tiers'] if t['minimum_pairs']==10)
    require({r['locus'] for r in pool}==set(expected),'Support-bin membership mismatch')
    old_context=rows(REVIEW/'gene_context.tsv')
    wanted={r['locus'] for r in pool}|{r['locus'] for r in old_context}
    catalog=PACKAGE/'references/mm10_te_fragments_2026-09-12/loci.tsv'
    loci={}
    with catalog.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['locus'] in wanted:
                r.update(start0=int(r['start0']),end=int(r['end']),chr=alias(r['chr']))
                loci[r['locus']]=r
    require(set(loci)==wanted,'Catalog missing loci')
    fai=Path(gb_path('references/mm10.fa.fai'))
    lengths={alias(v[0]):int(v[1]) for v in (line.split('\t') for line in fai.read_text().splitlines())}
    annotations={};inputs=[record(pq),record(gq),record(CHECK/'C6_pair_IAP_all.tsv'),record(CHECK/'C6K.support_scan.json'),record(catalog),record(fai),record(REVIEW/'gene_context.tsv')]
    for label,meta in gqc['annotations'].items():
        p=Path(meta['input']['path']);annotations[label]=parse(p,loci,lengths);inputs.append(record(p))
        for lid in {r['locus'] for r in old_context}:
            old=[r for r in old_context if r['locus']==lid and r['annotation']==label and int(r['gene_body_overlap_bp'])>0]
            actual=annotations[label][lid]['body']
            require(set(actual)=={r['gene_id'] for r in old},'Nine reference-gene memberships disagree')
            for r in old:
                g=actual[r['gene_id']]
                require(g['body']==int(r['gene_body_overlap_bp']) and g['exon']==int(r['exon_overlap_bp']),'Nine reference-gene overlap lengths disagree')
        print(now(),label,'Annotation and reference-gene checks passed',flush=True)
    output=[]
    for r in pool:
        lid=r['locus'];loc=loci[lid]
        out=dict(r);out['IGV_mm10']=f"{lid.split('@',1)[1].split(':',1)[0]}:{loc['start0']+1}-{loc['end']}"
        for label in ('Novogene','RefSeq'):
            ctx=annotations[label][lid]
            for field in ('category','body_genes','exon_genes','exon_bp'):out[label+'_'+field]=ctx[field]
        n=out['Novogene_category'];s=out['RefSeq_category']
        out['context_both_annotations']=n if n==s else 'annotation_dependent'
        m=float(r['C6M_final_count']);k=float(r['C6K_final_count'])
        out['KCl_higher_mapped_depth']=int(k*q['samples']['C6M']['mapped_pairs']>m*q['samples']['C6K']['mapped_pairs'])
        out['KCl_at_least_half_confident']=int(float(r['C6K_final_conf'])>=.5*k)
        out['KCl_support20']=int(float(r['C6K_final_conf'])>=20 and float(r['C6K_init_best'])>=20)
        output.append(out)
    output.sort(key=lambda r:(r['context_both_annotations'],-r['KCl_higher_mapped_depth'],-float(r['C6K_final_conf']),r['locus']))
    path=CHECK/(OUT+'.tsv');table(path,list(output[0]),[[r[c] for c in output[0]] for r in output])
    reread=rows(path)
    require(len(reread)==len(output) and {r['locus'] for r in reread}==set(expected),'Output membership mismatch')
    for a,b in zip(output,reread):require(all(('' if v is None else str(v))==b[k] for k,v in a.items()),'Output cell mismatch')
    subsets={'support10':output,
        'support10_KCl_higher':[r for r in output if r['KCl_higher_mapped_depth']],
        'support10_KCl_higher_half_confident':[r for r in output if r['KCl_higher_mapped_depth'] and r['KCl_at_least_half_confident']],
        'support20_KCl_higher':[r for r in output if r['KCl_higher_mapped_depth'] and r['KCl_support20']],
        'support20_KCl_higher_half_confident':[r for r in output if r['KCl_higher_mapped_depth'] and r['KCl_support20'] and r['KCl_at_least_half_confident']]}
    summary={name:summarize(rr) for name,rr in subsets.items()}
    require(all(sum(v['categories'].values())==v['n'] for v in summary.values()),'Non-disjoint classifications')
    data=dict(status='descriptive_context_complete',created_utc=now(),source=record(Path(__file__)),helper=record(PACKAGE/'R/telescope_common.py'),
        inputs=inputs,outputs=[record(path)],summary=summary,units='RepeatMasker fragments',
        QC=dict(source_hashes_verified=True,support_membership_verified=True,original_nine_gene_and_exon_overlaps_reproduced=True,output_cells_readback_verified=True,category_partition_verified=True),
        methods=['C6K support requires >=10 posterior-confident and >=10 sole-best-score fragments. Additional summaries use 20-fragment and half-confident bins.',
          'KCl-higher means counts per million total mapped fragments exceed C6M, without a pseudocount or additional effect-size or significance threshold.',
          'Uses the same Novogene GRCm38.p6 and ncbiRefSeq.2021-04-23 annotations; source coordinates converted from 1-based GTF to 0-based half-open BED.',
          'Gene bodies use gene/transcript envelopes. Fully_intronic requires complete containment in at least one body and zero exon overlap across all overlapping genes in that annotation.',
          'Outside_gene_bodies requires no body or exon overlap. Both annotations must agree for a shared class; disagreement stays explicit. Strand ignored because libraries are treated as unstranded.'])
    verify_records(inputs+data['outputs']+[data['source'],data['helper']])
    with (CHECK/(OUT+'.json')).open('x') as f:json.dump(data,f,indent=2);f.write('\n')
    require(json.loads((CHECK/(OUT+'.json')).read_text())==data,'Summary readback mismatch')
    print(json.dumps({k:{a:b for a,b in v.items() if a!='loci'} for k,v in summary.items()},indent=2))
if __name__=='__main__':main()
