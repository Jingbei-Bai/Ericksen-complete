"""Full source normal on the eight nonconic three-p RR3 candidates."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_single_RR3_normal.sobj'));OP=d['ring'];OK=OP.base_ring();K=PolynomialRing(QQ,names=('tn_k','tn_r')).fraction_field();k,r=K.gens();m=1+k;U=PolynomialRing(K,'tn_x');X=U.gen();F=X**8-2*m*X**6+2*m*m*X**4-m**3*X*X+4*m**4;E=U.quotient(F,'tn_b');x=E.gen();z=x-x**3/m;u=x*x/m;H=K(2);Z=K(-1);D=(H-Z)*(Z-1);e=x/Z*(-D*u**3+3*D*u**2+(-H*(Z-1)+3*Z*Z-Z)*u-2*(H*(Z-1)+Z*Z));images=(x,z,E.zero(),E.zero(),E.zero(),e,E.zero());ph=OK.ring().hom([k/r,r,H,Z],K);cache={}
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def mp(f):
 ans=E.zero()
 for ex,c in f.dict().items():
  ex=tuple(map(int,ex))
  if any(n and not v for n,v in zip(ex,images)):continue
  ans+=cv(c)*prod(v**n for n,v in zip(ex,images))
 return ans
M=matrix(E,28,9,[mp(f) for f in d['normal_matrix'].list()]);boundary=[mp(f) for f in d['boundary']];assert all(not f for f in boundary);print('BUILT',time.monotonic()-st,flush=True)
save({'ring':E,'matrix':M,'candidate_images':images,'boundary':boundary,'polynomial':F,'scope':__doc__},str(p/'independent_subcritical_single_RR_three_p_second_normal.sobj'),compress=False)

print('Fresh whole source matrix; unused generic rank diagnostic omitted',flush=True)
