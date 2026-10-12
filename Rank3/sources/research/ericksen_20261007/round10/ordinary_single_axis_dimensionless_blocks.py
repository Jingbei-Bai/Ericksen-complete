from sage.all import *
from pathlib import Path
import time,os,resource,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'ordinary_single_axis_J9.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
U=d['ring'];OK=U.base_ring();BK=PolynomialRing(QQ,names=('a','r','H','Z'));K=BK.fraction_field();a,r,H,Z=K.gens();q=1+r;m=1+a*r;P=PolynomialRing(K,names=('u','v'));u,v=P.gens();ph=OK.hom([a,a*q/m,H/r,q*Z/r],K);cache={}
def cv(c):
 if c not in cache:cache[c]=ph(c)
 return cache[c]
def mp(f):return P({tuple(map(int,e)):cv(c) for e,c in U(f).dict().items()})
A=d['matrix'].apply_map(mp);B=d['matrix_u_nonzero_seven'].apply_map(mp)
print('MAPPED',time.monotonic()-st,flush=True)
for i,row in enumerate(B.rows()):
 if any(row):print('SUPPORT',i,[j for j,c in enumerate(row) if c],flush=True)
# Linear range left elimination using only constant original source pivots.
C=matrix(P,B);ops=[];piv=[]
for c in range(5):
 candidates=[i for i in range(40) if i not in piv and C[i,c] and C[i,c].total_degree()==0]
 if not candidates:
  print('NOCONSTANT',c,flush=True);break
 ii=min(candidates,key=lambda i:sum(bool(x) for x in C.row(i)));pv=C[ii,c];piv.append(ii)
 print('PIVOT',ii,c,pv.factor(),flush=True)
 for i in range(40):
  if i==ii or not C[i,c]:continue
  cc=C[i,c]/pv;C.set_row(i,[P(C[i,j]-cc*C[ii,j]) for j in range(7)]);ops.append((i,ii,cc))
 # Normalize just chosen row, source pivot retained.
 C.set_row(ii,[P(x/pv) for x in C.row(ii)])
save({'ring':P,'matrix9':A,'matrix7':B,'reduced':C,'pivot_rows':piv,'operations':ops,'beta':mp(d['beta']),'lambda':mp(d['lambda'])},str(p/'ordinary_single_axis_dimensionless_blocks.sobj'),compress=False)
print('REDUCED',len(piv),time.monotonic()-st,flush=True)
for i,row in enumerate(C.rows()):
 if i not in piv and any(row):print('REMAIN',i,[(j,str(x.factor())) for j,x in enumerate(row) if x],flush=True)
