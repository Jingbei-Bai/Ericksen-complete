"""Exact generic-base quadratic range correction on m*l^2=-ar.
Whole coupled columns are retained; the cyclic correction is not extended
by assumption. Any nonzero residual is saved as an actual obstruction.
"""
from sage.all import *
from pathlib import Path
import json,time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'ordinary_pair_double_kernel_coupled.sobj'));P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();q=1+r;k=a*r;m=1+k;l=P.gen();M=d['matrix'];labs=d['labels'];mu=m*l*l+k
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
def ck(n,x):assert bool(x),n;checks.append(n)
def rem(f):return sum((c*(-k/m)**(int(e)//2)*l**(int(e)%2) for e,c in P(f).dict().items()),P.zero())
def inv(f):
    f=rem(f);bar=f(-l);n=rem(f*bar);ck('quadratic inverse norm scalar',n.degree()<=0 and n!=0)
    return P(bar/n[0])
ix=[labs.index(('old',i)) for i in (2,17)];J=M.matrix_from_rows_and_columns(ix,(1,2));Q=vector(P,[M[i,8] for i in ix]);D=rem(J.det());Di=inv(D)
X=[rem((-J[1,1]*Q[0]+J[0,1]*Q[1])*Di),rem((J[1,0]*Q[0]-J[0,0]*Q[1])*Di)]
for i,x in enumerate(X):print('RANGE',i,'factor',x.factor(),flush=True)
res=[rem(M[i,1]*X[0]+M[i,2]*X[1]+M[i,8]) for i in range(M.nrows())]
for i,z in enumerate(res):
    if z:print('RESIDUAL',labs[i],str(z.factor()),flush=True)
print('NONZERO residuals',sum(bool(z) for z in res),'secs',time.monotonic()-st,flush=True)
save({'ring':P,'pivot_rows':(2,17),'pivot_columns':(1,2),'pivot_matrix':J,'pivot_det_mod_wall':D,'correction':X,'residuals':res,'labels':labs,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_minus_general_quadratic_range.sobj'),compress=False)
(p/'ordinary_pair_minus_general_quadratic_range.json').write_text(json.dumps({'nonzero_residuals':sum(bool(z) for z in res),'seconds':time.monotonic()-st,'checks':len(checks)},indent=2)+'\n')
