"""Decode the shared B1 HMM and summarize LAD segments and replicate agreement."""
from pathlib import Path
import argparse
import itertools
import json
import numpy as np
import pandas as pd

GROUPS=[[0,1],[2,3,4],[5,6,7]]
NAMES=['TTX','Bicuculline','Bicuculline + nimodipine']

def close(x,y,atol=1e-10):
    np.testing.assert_allclose(x,y,atol=atol,rtol=1e-10)

def decode(model,x,chrom,start):
    """Gaussian Viterbi decoding, reset at chromosomes and support gaps."""
    emit=[j for j,s in enumerate(model['states']) if s['distribution'] is not None]
    mu=np.array([model['states'][j]['distribution']['parameters'][0] for j in emit])
    sd=np.array([model['states'][j]['distribution']['parameters'][1] for j in emit])
    rank=np.argsort(np.argsort(mu,kind='stable'),kind='stable')
    tr=np.zeros((3,3));initial=np.zeros(3)
    for a,b,p,*unused in model['edges']:
        if a==model['start_index'] and b in emit:initial[emit.index(b)]=p
        if a in emit and b in emit:tr[emit.index(a),emit.index(b)]=p
    close(tr.sum(1),1);close(initial.sum(),1)
    cuts=np.r_[0,np.flatnonzero((chrom[1:]!=chrom[:-1])|(start[1:]!=start[:-1]+30000))+1,len(x)]
    out=np.empty(x.shape,dtype=np.int8)
    with np.errstate(divide='ignore'):lt=np.log(tr);li=np.log(initial)
    # Decode all libraries together, resetting at chromosome and support gaps.
    for lo,hi in zip(cuts[:-1],cuts[1:]):
        vals=x[lo:hi];score=li[None,:]-np.log(sd*np.sqrt(2*np.pi))-.5*((vals[0,:,None]-mu)/sd)**2
        back=np.zeros((hi-lo,8,3),dtype=np.int8)
        for t in range(1,hi-lo):
            z=score[:,:,None]+lt[None,:,:];back[t]=np.argmax(z,axis=1)
            score=np.max(z,axis=1)-np.log(sd*np.sqrt(2*np.pi))-.5*((vals[t,:,None]-mu)/sd)**2
        state=np.argmax(score,axis=1);out[hi-1]=rank[state]
        for t in range(hi-lo-1,0,-1):
            state=back[t,np.arange(8),state];out[lo+t-1]=rank[state]
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path)
    parser.add_argument('--hmm-input',type=Path)
    parser.add_argument('--signals',type=Path)
    parser.add_argument('--model',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();source=args.inputs;output=args.output
    output.mkdir(parents=True,exist_ok=True)
    hmm_input=args.hmm_input or (source/'shared_hmm_input.npz' if source else None)
    signals=args.signals or (source/'b1_signal_vectors.npz' if source else None)
    model_path=args.model or (source/'selected_model.json' if source else None)
    if not all((hmm_input,signals,model_path)):
        parser.error('Provide --inputs or explicit --hmm-input, --signals and --model paths.')
    frozen=dict(np.load(hmm_input))
    signal=dict(np.load(signals))
    model=json.loads(model_path.read_text())
    states=decode(model,frozen['x'],frozen['chromosome'],frozen['start'])
    np.savez_compressed(output/'decoded_states.npz',states=states)
    bp=signal['clean_bp'];chrom=signal['chromosome'];start=signal['start']
    state=np.column_stack([states[:,group].max(1) for group in GROUPS])
    signature=((state>=1)*[1,2,4]).sum(1)
    rows=[]
    for i,j in itertools.product(range(8),repeat=2):
        a=states[:,i]>=1;b=states[:,j]>=1
        rows.append(dict(sample1=signal['IP_ids'][i],sample2=signal['IP_ids'][j],Jaccard=bp[a&b].sum()/bp[a|b].sum()))
    pd.DataFrame(rows).to_csv(output/'replicate_Jaccard.tsv',sep='\t',index=False)
    sharing=[dict(signature=s,condition=name,state=k,clean_bp=int(bp[(signature==s)&(state[:,j]==k)].sum()))
             for s in range(8) for j,name in enumerate(NAMES) for k in range(3)]
    pd.DataFrame(sharing).to_csv(output/'nested_sharing.tsv',sep='\t',index=False)
    members=[];owner=np.full(len(signature),-1,int)
    for i in np.flatnonzero(signature):
        if members and i==members[-1][-1]+1 and chrom[i]==chrom[i-1] and start[i]==start[i-1]+30000 and signature[i]==signature[i-1]:
            members[-1].append(i)
        else:
            members.append([i])
    regions=[]
    for k,indices in enumerate(members):
        owner[indices]=k
        regions.append(dict(region_index=k,chromosome=chrom[indices[0]],start0=start[indices[0]],
                            end0=start[indices[-1]]+30000,clean_bp=bp[indices].sum(),
                            signature=signature[indices[0]],bins=len(indices)))
    values=np.array([np.average(signal['condition_depth'][indices],axis=0,weights=bp[indices]) for indices in members])
    matrix=values.T-values.mean(1);u,d,vh=np.linalg.svd(matrix,full_matrices=False);scores=u*d
    for j in range(scores.shape[1]):
        scores[:,j]*=np.sign(vh[j,np.argmax(abs(vh[j]))])
    pd.DataFrame(scores,columns=[f'PC{i+1}' for i in range(8)]).assign(sample=signal['IP_ids']).to_csv(output/'domain_PCA_scores.tsv',sep='\t',index=False)
    pd.DataFrame(dict(component=np.arange(1,9),percent=d*d/(d*d).sum()*100)).to_csv(output/'domain_PCA_variance.tsv',sep='\t',index=False)
    pd.DataFrame(regions).to_csv(output/'all_condition_regions.tsv',sep='\t',index=False)
    np.savez_compressed(output/'all_condition_regions.npz',owner=owner,signal=values,signature=signature,states=state)
    print(f'Wrote {len(regions):,} LAD segments and summaries for {states.shape[1]} B1 libraries.')

if __name__=='__main__':
    main()
