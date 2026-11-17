"""Count eligible paired fragments in 5-kb bins on chrX and chrY."""
from pathlib import Path
import argparse
import bisect
import csv
import subprocess
import numpy as np
import pysam

SEX=['chrX','chrY']

def table(path, rows):
    with path.open('w', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader();writer.writerows(rows)

def overlap(a,b,start,end):
    i=bisect.bisect_left(a,end)-1
    return i>=0 and b[i]>start

def prepare_reference(fasta, mask, output):
    intervals={c:[] for c in SEX}
    for line in mask.read_text().splitlines():
        if line.startswith(('#','track','browser')) or not line:continue
        c,a,b=line.split()[:3]
        if c in intervals:intervals[c].append((int(a),int(b)))
    out={};rr=[]; lengths=[]
    with pysam.FastaFile(str(fasta)) as f:
        for c in SEX:
            seq=np.frombuffer(f.fetch(c).upper().encode(),dtype=np.uint8);length=len(seq);lengths.append(length)
            acgt=np.isin(seq,[65,67,71,84]);good=acgt.copy();merged=[]
            assert intervals[c],c+' mask unexpectedly absent'
            for a,b in sorted(intervals[c]):
                assert 0<=a<b<=length
                if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
                else:merged.append((a,b))
            for a,b in merged:good[a:b]=False
            ss=np.array([a for a,b in merged]);ee=np.array([b for a,b in merged]);st=np.arange(0,length,5000)
            clean=np.add.reduceat(good,st).astype(np.int64)
            change=np.diff(np.r_[False,~acgt,False].astype(np.int8))
            out.update({c+'_clean_bp':clean,c+'_actual_bp':np.minimum(5000,length-st),
                c+'_mask_start':ss,c+'_mask_end':ee,c+'_nonACGT_start':np.flatnonzero(change==1),c+'_nonACGT_end':np.flatnonzero(change==-1)})
            for k in np.unique(np.r_[np.arange(0,len(st),101),len(st)-1]):
                lo=int(st[k]);hi=min(lo+5000,length)
                independent=int(acgt[lo:hi].sum())
                for a,b in merged:
                    if a>=hi:break
                    if b>lo:independent-=int(acgt[max(a,lo):min(b,hi)].sum())
                assert independent==clean[k]
            rr.append({'chromosome':c,'length_bp':length,'clean_bp':int(clean.sum()),'mask_intervals':len(merged)})
    out['chromosomes']=np.array(SEX);out['lengths']=np.array(lengths)
    np.savez_compressed(output/'XY_reference.npz',**out);table(output/'XY_reference_summary.tsv',rr)

def count_sample(sid, bam, reference, output):
    subprocess.run(['samtools','quickcheck','-v',str(bam)],check=True)
    ref=np.load(reference);out={};rr=[]
    with pysam.AlignmentFile(str(bam),'rb',threads=1) as f:
        for c,length in zip(ref['chromosomes'],ref['lengths']):
            c=str(c);length=int(length);clean_bp=ref[c+'_clean_bp'];ct=np.zeros(len(clean_bp),dtype=np.int64);mate=ct.copy()
            ss,ee=ref[c+'_mask_start'].tolist(),ref[c+'_mask_end'].tolist()
            ns,ne=ref[c+'_nonACGT_start'].tolist(),ref[c+'_nonACGT_end'].tolist()
            eligible=masked=nonacgt=0
            for read in f.fetch(c):
                flag=read.flag
                if not flag&2 or flag&(4|8|256|512|1024|2048):continue
                if not (flag&64 or flag&128):continue
                if read.next_reference_id!=read.reference_id or read.template_length==0:continue
                assert read.has_tag('MQ')
                if min(read.mapping_quality,read.get_tag('MQ'))<10:continue
                read1=bool(flag&64);a=min(read.reference_start,read.next_reference_start);b=a+abs(read.template_length)
                assert 0<=a<b<=length
                mid=a+(b-a)//2;k=mid//5000
                if read1:eligible+=1
                if overlap(ss,ee,a,b):
                    if read1:masked+=1
                    continue
                if overlap(ns,ne,mid,mid+1):
                    if read1:nonacgt+=1
                    continue
                assert clean_bp[k]>0
                if read1:ct[k]+=1
                else:mate[k]+=1
            assert np.array_equal(ct,mate),(sid,c,'independent mate counts')
            assert int(ct.sum())==eligible-masked-nonacgt
            out[c]=ct
            rr.append(dict(sample=sid,chromosome=c,eligible_fragments=eligible,mask_excluded=masked,
                nonACGT_midpoints_excluded=nonacgt,retained=int(ct.sum()),mate_recount_exact=True))
    np.savez_compressed(output/(sid+'.XY_counts.npz'),**out)
    return rr

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    ref=sub.add_parser('reference')
    ref.add_argument('--fasta',type=Path,required=True)
    ref.add_argument('--mask',type=Path,required=True)
    ref.add_argument('--output',type=Path,required=True)
    count=sub.add_parser('count')
    count.add_argument('--sample',required=True)
    count.add_argument('--bam',type=Path,required=True)
    count.add_argument('--reference',type=Path,required=True)
    count.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if args.command=='reference':
        prepare_reference(args.fasta,args.mask,args.output)
    else:
        rows=count_sample(args.sample,args.bam,args.reference,args.output)
        table(args.output/(args.sample+'.XY_counts.tsv'),rows)

if __name__=='__main__':
    main()
