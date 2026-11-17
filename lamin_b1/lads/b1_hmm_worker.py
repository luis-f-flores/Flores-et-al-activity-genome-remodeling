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

import os
os.environ.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
import sys, csv, gzip, json, random, time, itertools
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
import pomegranate
from pomegranate import HiddenMarkovModel, NormalDistribution
assert pomegranate.__version__ == '0.11.1'

def segments(ch, start, end):
    assert len(ch) and len(ch) == len(start) == len(end)
    assert np.all(end > start)
    cuts = np.r_[0, np.flatnonzero((ch[1:] != ch[:-1]) | (start[1:] != end[:-1]))+1, len(ch)]
    return [(int(a), int(b)) for a,b in zip(cuts[:-1],cuts[1:])]

def parameters(model):
    emitting = [i for i,s in enumerate(model.states) if s.distribution is not None]
    assert emitting == list(range(len(emitting))), 'Check posterior-column mapping for this runtime'
    mu = np.array([model.states[i].distribution.parameters[0] for i in emitting])
    sd = np.array([model.states[i].distribution.parameters[1] for i in emitting])
    assert np.isfinite(mu).all() and np.isfinite(sd).all() and (sd>0).all()
    order = np.argsort(mu, kind='stable')
    assert (np.diff(mu[order]) > 0).all(), 'Indistinguishable state means'
    rank = {emitting[j]:i for i,j in enumerate(order)}
    mat = model.dense_transition_matrix()
    initial = mat[model.states.index(model.start), emitting]
    trans = mat[np.ix_(emitting,emitting)]
    assert np.allclose(initial.sum(),1) and np.allclose(trans.sum(1),1)
    return emitting, mu, sd, order, rank, initial, trans

def fit(sequences, k):
    assert all(np.ndim(v)==1 and len(v) and np.isfinite(v).all() for v in sequences)
    np.random.seed(12345); random.seed(12345)
    return HiddenMarkovModel.from_samples(NormalDistribution, n_components=k,
        X=sequences, algorithm='baum-welch', n_jobs=1, n_init=1, init='kmeans++',
        max_kmeans_iterations=1, max_iterations=1000, stop_threshold=1e-9,
        min_iterations=0, end_state=False, return_history=True, verbose=False)

def decode(model, sequences):
    emitting,mu,sd,order,rank,initial,trans = parameters(model)
    paths=[]; probs=[]
    for x in sequences:
        raw=model.predict(x, algorithm='viterbi')
        path=np.array([rank[i] for i in raw if i in rank],dtype=np.int8)
        p=model.predict_proba(x)[:,order]
        assert len(path)==len(x) and p.shape==(len(x),len(emitting))
        assert np.isfinite(p).all() and np.allclose(p.sum(1),1,atol=1e-7)
        paths.append(path); probs.append(p)
    return np.concatenate(paths),np.concatenate(probs)

def independent(model,x):
    emitting,mu,sd,order,rank,initial,trans=parameters(model)
    emit=-np.log(sd*np.sqrt(2*np.pi))-.5*((x[:,None]-mu)/sd)**2
    with np.errstate(divide='ignore'): li=np.log(initial); lt=np.log(trans)
    a=np.zeros_like(emit); b=np.zeros_like(emit); v=np.zeros_like(emit); back=np.zeros_like(emit,dtype=int)
    a[0]=li+emit[0]; v[0]=a[0]
    for i in range(1,len(x)):
        a[i]=emit[i]+logsumexp(a[i-1,:,None]+lt,axis=0)
        cand=v[i-1,:,None]+lt
        back[i]=np.argmax(cand,axis=0);v[i]=emit[i]+np.max(cand,axis=0)
    for i in range(len(x)-2,-1,-1): b[i]=logsumexp(lt+emit[i+1]+b[i+1],axis=1)
    lp=float(logsumexp(a[-1])); post=np.exp(a+b-lp)
    path=[int(np.argmax(v[-1]))]
    for i in range(len(x)-1,0,-1): path.append(int(back[i,path[-1]]))
    path=np.array([rank[i] for i in path[::-1]])
    actual,probs=decode(model,[x])
    err=float(np.max(abs(probs-post[:,order])))
    assert err<1e-7 and np.array_equal(path,actual)
    assert abs(lp-model.log_probability(x)) < 1e-7*max(1,abs(lp))
    return {'posterior_max_error':err,'Viterbi_exact':True,'loglik_replay':True}

def domains(ch,st,en,state,target):
    out=[]
    for i in np.flatnonzero(state==target):
        if out and out[-1][0]==ch[i] and out[-1][2]==int(st[i]):out[-1][2]=int(en[i])
        else:out.append([str(ch[i]),int(st[i]),int(en[i])])
    return out

def self_test():
    records=[]
    
    ch=np.array(['chr1']*3+['chrX']);st=np.array([0,20,60,0]);en=st+20
    assert segments(ch,st,en)==[(0,2),(2,3),(3,4)]
    assert domains(ch,st,en,np.ones(4),1)==[['chr1',0,40],['chr1',60,80],['chrX',0,20]]
    for k in [2,3]:
        rng=np.random.RandomState(130+k)
        seq=[np.concatenate([rng.normal(1+3*j,.25,80) for j in range(k)]) for _ in range(3)]
        truth=np.tile(np.repeat(np.arange(k),80),3)
        m,h=fit(seq,k);path,probs=decode(m,seq)
        assert np.mean(path==truth)>.99
        repeat,_=fit(seq,k);assert np.array_equal(path,decode(repeat,seq)[0])
        restored=HiddenMarkovModel.from_json(m.to_json())
        assert np.allclose(probs,decode(restored,seq)[1],atol=1e-12)
        checks=[independent(m,x) for x in [seq[0][:200],np.array([2.5]),np.array([1.2,2.5,4.2,1.2])]]
        
        e,mu,sd,order,rank,ini,tr=parameters(m); x=np.array([1.1,2.3,4.1,1.4])
        ll=[]
        with np.errstate(divide='ignore'):
            for states in itertools.product(range(k),repeat=len(x)):
                score=np.log(ini[states[0]])
                for i,s in enumerate(states):
                    score+=-np.log(sd[s]*np.sqrt(2*np.pi))-.5*((x[i]-mu[s])/sd[s])**2
                    if i:score+=np.log(tr[states[i-1],s])
                ll.append(score)
        assert abs(logsumexp(ll)-m.log_probability(x))<1e-8
        
        solo=decode(m,[seq[0]])
        combined=decode(m,[seq[0],np.array([1.,4.,7.])])
        assert np.array_equal(solo[0],combined[0][:len(seq[0])])
        assert np.array_equal(solo[1],combined[1][:len(seq[0])])
        records.append({'states':k,'checks':checks,'seed_repeat':True,'roundtrip':True,
            'brute_force_loglik':True,'XY_does_not_change_autosomal_decoding':True})
    print(json.dumps({'status':'PASS','coordinate_and_gap_checks':True,'checks':records},indent=2))

def main():
    source,sid,prefix,k=sys.argv[1:];k=int(k);starttime=time.monotonic()
    d=np.genfromtxt(source,delimiter='\t',names=True,dtype=None,encoding='utf-8')
    ch=d['chr'];st=d['start'];en=d['end'];x=d[sid].astype(float)
    assert np.isfinite(x).all() and (x>0).all()
    seqs=segments(ch,st,en)
    training=[x[a:b] for a,b in seqs if str(ch[a]) in ['chr'+str(j) for j in range(1,20)]]
    model,history=fit(training,k)
    emit,mu,sd,order,rank,init,trans=parameters(model)
    state,probs=decode(model,[x[a:b] for a,b in seqs])
    assert len(state)==len(x)
    checks=[independent(model,x[a:min(b,a+200)]) for a,b in [seqs[0],max(seqs,key=lambda s:s[1]-s[0]),seqs[-1]]]
    restored=HiddenMarkovModel.from_json(model.to_json())
    a,b=max(seqs,key=lambda s:s[1]-s[0]);test=x[a:min(b,a+200)]
    assert np.array_equal(decode(model,[test])[0],decode(restored,[test])[0])
    assert np.allclose(decode(model,[test])[1],decode(restored,[test])[1],atol=1e-12)
    with gzip.open(prefix+'.posterior.tsv.gz','wt') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['row_id','state','p_LAD']+['p_state'+str(j) for j in range(k)])
        w.writerows([i,int(state[i]),probs[i,-1]]+probs[i].tolist() for i in range(len(x)))
    names=['low','high'] if k==2 else ['low','intermediate','high']
    for j,name in enumerate(names):
        with open(prefix+'.'+name+'.bed','w') as f:
            csv.writer(f,delimiter='\t').writerows(domains(ch,st,en,state,j))
    Path(prefix+'.model.json').write_text(model.to_json())
    ll=np.array(history.log_probabilities);delta=np.array(history.improvements)
    assert np.isfinite(ll).all() and np.isfinite(delta).all()
    tol=max(1e-7,1e-10*max(1,abs(float(ll[-1]))))
    converged=bool(len(ll)<1000 and delta[-1]<1e-9 and delta.min()>=-tol)
    hist=[{'iteration':int(i),'loglik':float(l),'improvement':float(v)} for i,l,v in zip(history.epochs,ll,delta)]
    with open(prefix+'.iterations.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(hist[0]),delimiter='\t');w.writeheader();w.writerows(hist)
    met={'sample':sid,'states':k,'engine':'pomegranate','smoothing':'N/A','bins':len(x),
        'training_bins':sum(map(len,training)),'training_sequences':len(training),'prediction_sequences':len(seqs),
        'iterations':len(ll),'converged':converged,'last_improvement':float(delta[-1]),'min_improvement':float(delta.min()),
        'means_ordered':mu[order].tolist(),'SD_ordered':sd[order].tolist(),
        'state_bins':np.bincount(state,minlength=k).tolist(),'seconds':time.monotonic()-starttime,
        'seed':12345,'input':source,'QC':checks,'model_roundtrip':True}
    Path(prefix+'.metrics.json').write_text(json.dumps(met,indent=2))
    print(json.dumps(met),flush=True)

if __name__=='__main__':
    if sys.argv[1:]==['--self-test']:self_test()
    else:main()
