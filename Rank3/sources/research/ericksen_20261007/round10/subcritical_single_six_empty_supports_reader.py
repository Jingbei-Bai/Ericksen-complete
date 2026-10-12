"""Actual original-source row certificates for the six empty single-RR supports.

The cyclic row identities are recomputed from original full boundary rows;
we do not specialize a generic echelon transformation with unknown pivots.
"""
from sage.all import *
from pathlib import Path
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_inactive_single_RR_boundary.sobj'));OP=d['ring'];OK=OP.base_ring();K=PolynomialRing(QQ,names=('k','r')).fraction_field();k,r=K.gens();P=PolynomialRing(K,names=OP.variable_names());vs=P.gens();x,y,z,w,e0,e1,e2,de=vs;ph=OK.ring().hom([k/r,r,1,-1],K)
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2)
def mp(f):return P({tuple(map(int,ee)):ph(c.numerator())/ph(c.denominator()) for ee,c in f.dict().items()})
src=list(map(mp,d['items'][2]['rows']));ev=list(map(mp,d['items'][2]['even']));checks=[];items=[]
def ck(n,v):assert bool(v),n;checks.append(n)
def protected(c):
 den=K(c).denominator()
 for f,n in den.factor():
  if f.total_degree()==0:continue
  lc=f.leading_coefficient();ff=f/lc
  if ff not in (k,r,k+1,r+1):return False
 return True
targets={3:e0*y,6:e0*y,9:e2*x,12:e2*w,11:x*y*w,14:y*z*w}
for mask,target in targets.items():
 support=[j for j in range(4) if mask>>j&1];sm=P.hom([vs[j] if j in support else P.zero() for j in range(4)]+[e0,P.zero(),e2,P.zero()],P);qq=sm(src[6]);ck('Q0 e1 independent unit mask '+str(mask),src[6].derivative(e1).constant_coefficient()!=0);ck('Q0 zero after e1 elimination '+str(mask),qq==0)
 candidates=[];orig=[]
 for i,f in enumerate(src):
  f=sm(f)
  if not f:continue
  for j in support:
   if all(ee[j] for ee in f.dict()):
    nf=P({tuple(int(n)-(1 if a==j else 0) for a,n in enumerate(ee)):c for ee,c in f.dict().items()});candidates.append(nf);orig.append((i,j));ck('literal active-p division '+str((mask,i,j)),nf*vs[j]==f)
 mons=sorted(set(target.dict()).union(*(set(f.dict()) for f in candidates)));A=matrix(K,[[f.monomial_coefficient(P({ee:1})) for f in candidates] for ee in mons]);tv=vector(K,[target.monomial_coefficient(P({ee:1})) for ee in mons]);co=A.solve_right(tv)
 ck('actual original cyclic row combination '+str(mask),sum((c*f for c,f in zip(co,candidates)),P.zero())==target)
 ck('all material cofactor denominators protected '+str(mask),all(protected(c) for c in co))
 ee=list(map(sm,ev));res=None
 if mask==3:ck('mask3 actual Q6 e2 forcing',ee[2]==-k/(2*(1+k))*x*y*e2)
 if mask==6:ck('mask6 actual Q6 e2 forcing',ee[2]==k/(2*(1+k))*z*y*(10*y*z+e2))
 if mask==9:ck('mask9 actual Q6 e0 forcing',ee[2]==1/(2*(1+k))*w*x*(10*x*w+e0))
 if mask==12:ck('mask12 actual Q6 e0 forcing',ee[2]==-1/(2*(1+k))*e0*w*z)
 if mask==3:res=P.hom([x,y,z,w,P.zero(),P.zero(),P.zero(),P.zero()],P)(ee[0]);ck('mask3 Q6 and Q1 contradiction',res==8*x*x*y*y)
 if mask==6:res=P.hom([x,y,z,w,P.zero(),P.zero(),-10*y*z,P.zero()],P)(ee[0]);ck('mask6 Q6 and Q1 contradiction',res==16*y*y*z*z)
 if mask==9:res=P.hom([x,y,z,w,-10*x*w,P.zero(),P.zero(),P.zero()],P)(ee[1]);ck('mask9 Q6 and Q5 contradiction',res==16*x*x*w*w)
 if mask==12:res=P.hom([x,y,z,w,P.zero(),P.zero(),P.zero(),P.zero()],P)(ee[1]);ck('mask12 Q6 and Q5 contradiction',res==8*z*z*w*w)
 items.append({'mask':mask,'source_references':orig,'cofactors':co,'target':target,'even_residual':res});print('MASK',mask,'PASS',flush=True)
save({'ring':P,'items':items,'checks':checks,'scope':__doc__},str(p/'subcritical_single_six_empty_supports_reader.sobj'),compress=False);(p/'subcritical_single_six_empty_supports_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS',len(checks),'sec',time.monotonic()-st,flush=True)
