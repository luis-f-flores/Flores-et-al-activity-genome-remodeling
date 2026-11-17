"""LAD-body and oriented-boundary query windows."""
import numpy as np

def profile_queries(regions,kind):
    q=[]
    if kind=='body':
        xpos=np.r_[(np.linspace(-.5,0,18)[:-1]+np.linspace(-.5,0,18)[1:])/2,(np.linspace(0,1,51)[:-1]+np.linspace(0,1,51)[1:])/2,(np.linspace(1,1.5,18)[:-1]+np.linspace(1,1.5,18)[1:])/2]
        for rid,(c,a,b) in enumerate(regions):
            edges=np.r_[np.linspace(a-500000,a,18),np.linspace(a,b,51)[1:],np.linspace(b,b+500000,18)[1:]]
            for pos,(lo,hi) in enumerate(zip(edges[:-1],edges[1:])):q.append((c,float(lo),float(hi),pos,rid))
    else:
        edges=np.linspace(-500000,500000,35);xpos=(edges[:-1]+edges[1:])/2000
        for rid,(c,a,b) in enumerate(regions):
            for side in [0,1]:
                for pos,(x0,x1) in enumerate(zip(edges[:-1],edges[1:])):
                    if side==0:lo=a+x0;hi=min(a+x1,b) if x1>0 else a+x1
                    else:lo=max(b-x1,a) if x1>0 else b-x1;hi=b-x0
                    q.append((c,float(lo),float(hi),pos,2*rid+side))
    return q,xpos
