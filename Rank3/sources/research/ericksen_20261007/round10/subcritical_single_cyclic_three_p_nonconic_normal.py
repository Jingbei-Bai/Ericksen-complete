"""Full source normal on the eight nonconic three-p RR3 candidates."""
from sage.all import *
from pathlib import Path
import time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_single_RR3_normal.sobj'));OP=d['ring'];OK=OP.base_ring();K=PolynomialRing(QQ,names=('tn_k','tn_r')).fraction_field();k,r=K.gens();m=1+k;U=PolynomialRing(K,'tn_x');X=U.gen();F=X**8-2*m*X**6+2*m*m*X**4-m**3*X*X+4*m**4;E=U.quotient(F,'tn_b');x=E.gen();z=x-x**3/m;u=x*x/m;H=K(1);Z=K(-1);D=(H-Z)*(Z-1);e=x/Z*(-D*u**3+3*D*u**2+(-H*(Z-1)+3*Z*Z-Z)*u-2*(H*(Z-1)+Z*Z));images=(x,z,E.zero(),E.zero(),E.zero(),e,E.zero());ph=OK.ring().hom([k/r,r,H,Z],K);cache={}
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2)

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
save({'ring':E,'matrix':M,'candidate_images':images,'boundary':boundary,'polynomial':F,'scope':__doc__},str(p/'subcritical_single_cyclic_three_p_nonconic_normal.sobj'),compress=False)

J=M.matrix_from_rows_and_columns(range(25),(2,3,4,6,8))
print('J5 RANK',J.rank(),'KERNEL',J.right_kernel_matrix(),flush=True)
save({'ring':E,'matrix':M,'J5matrix':J,'candidate_images':images,'polynomial':F},str(p/'subcritical_single_cyclic_three_p_nonconic_normal.sobj'),compress=False)
