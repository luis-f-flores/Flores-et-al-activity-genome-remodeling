"""Shared three-state B1 HMM training with balanced condition weighting."""
from pathlib import Path
import argparse,csv,json,random,time,sys,importlib.util
import numpy as np
MULT=[3,3,2,2,2,2,2,2]
def worker():
    spec=importlib.util.spec_from_file_location('gb_worker',Path(__file__).with_name('b1_hmm_worker.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def table(p,rows):
    with Path(p).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def fit_one(seed):
    m=worker();import pomegranate
    a=dict(np.load(R/'tables/frozen_input.npz'));x=a['x'];segments=m.segments(a['chromosome'],a['start'],a['end'])
    train=[(l,u) for l,u in segments if a['chromosome'][l] in ['chr'+str(i) for i in range(1,20)]]
    seq=[x[l:u,j] for j in range(8) for rep in range(MULT[j]) for l,u in train]
    np.random.seed(seed);random.seed(seed);t0=time.monotonic()
    print('FIT START seed',seed,'sequences',len(seq),'weighted bins',sum(map(len,seq)),flush=True)
    model,hist=m.HiddenMarkovModel.from_samples(m.NormalDistribution,n_components=3,X=seq,
      algorithm='baum-welch',n_jobs=1,n_init=1,init='kmeans++',max_kmeans_iterations=1,
      max_iterations=1000,stop_threshold=1e-9,min_iterations=0,end_state=False,return_history=True,verbose=False)
    prefix=R/'models'/('seed_'+str(seed));Path(str(prefix)+'.model.json').write_text(model.to_json())
    ll=np.asarray(hist.log_probabilities);delta=np.asarray(hist.improvements)
    assert np.isfinite(ll).all() and np.isfinite(delta).all()
    tol=max(1e-7,1e-10*max(1,abs(float(ll[-1]))))
    conv=bool(len(ll)<1000 and delta[-1]<1e-9 and delta.min()>=-tol)
    table(str(prefix)+'.iterations.tsv',[dict(iteration=int(i),loglik=float(v),improvement=float(d)) for i,v,d in zip(hist.epochs,ll,delta)])
    states=[];probs=[];checks=[]
    for j in range(8):
        st,pr=m.decode(model,[x[l:u,j] for l,u in segments]);states.append(st);probs.append(pr)
        l,u=max(segments,key=lambda z:z[1]-z[0]);checks.append(m.independent(model,x[l:min(u,l+200),j]))
    states=np.column_stack(states);probs=np.stack(probs,axis=1)
    restored=m.HiddenMarkovModel.from_json(model.to_json())
    l,u=max(segments,key=lambda z:z[1]-z[0]);test=x[l:min(u,l+300),0]
    np.testing.assert_array_equal(m.decode(model,[test])[0],m.decode(restored,[test])[0])
    np.testing.assert_allclose(m.decode(model,[test])[1],m.decode(restored,[test])[1],atol=1e-12,rtol=0)
    e,mu,sd,o,rank,init,tr=m.parameters(model)
    final_ll=sum(MULT[j]*sum(model.log_probability(x[l:u,j]) for l,u in train) for j in range(8))
    np.savez_compressed(str(prefix)+'.decoded.npz',states=states,posteriors=probs)
    metrics=dict(seed=seed,engine=pomegranate.__version__,python=sys.version,iterations=len(ll),converged=conv,
      last_improvement=float(delta[-1]),min_improvement=float(delta.min()),final_balanced_loglik=float(final_ll),
      last_history_loglik=float(ll[-1]),means_ordered=mu[o].tolist(),SD_ordered=sd[o].tolist(),
      transition_ordered=tr[np.ix_(o,o)].tolist(),initial_ordered=init[o].tolist(),
      state_bins_by_library=[np.bincount(states[:,j],minlength=3).tolist() for j in range(8)],
      independent_dynamic_programming=checks,model_roundtrip=True,seconds=time.monotonic()-t0)
    dump(str(prefix)+'.metrics.json',metrics);print(json.dumps(metrics),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--seeds',nargs='+',type=int,default=[12345,23456,34567]);a=p.parse_args();R=a.output
    assert not R.exists(),'Choose a new output directory.'
    (R/'models').mkdir(parents=True);(R/'tables').mkdir()
    import shutil
    shutil.copyfile(a.inputs,R/'tables/frozen_input.npz')
    for seed in a.seeds:fit_one(seed)
    metrics=[json.loads((R/'models'/('seed_'+str(seed)+'.metrics.json')).read_text()) for seed in a.seeds]
    valid=[m for m in metrics if m['converged']];assert valid,'No converged fit.'
    best=max(valid,key=lambda m:m['final_balanced_loglik']);print('Selected seed:',best['seed'])
