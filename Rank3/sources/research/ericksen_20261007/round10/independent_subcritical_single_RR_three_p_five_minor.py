"""True five-source minor, avoiding algebraic-field Gaussian inversion."""
from sage.all import *
from pathlib import Path
from itertools import permutations
import time,sys
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_three_p_five_module.sobj'));E=d['ring'];K=E.base_ring();k,r=K.gens();M=d['matrix'];rows=tuple(map(int,sys.argv[1].split(','))) if len(sys.argv)>1 else (2,3,6,21,9)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 n,den=c.numerator(),c.denominator();gg=n.gcd(den);n,den=n//gg,den//gg;lc=den.leading_coefficient();return K(n/lc)/K(den/lc)
def ne(v):return E([nk(v[i]) for i in range(8)])
A=M.matrix_from_rows(rows);det=E.zero();terms=0
for perm in permutations(range(5)):
 if any(not A[i,perm[i]] for i in range(5)):continue
 term=E.one();sign=(-1)**sum(perm[i]>perm[j] for i in range(5) for j in range(i+1,5))
 for i in range(5):term=ne(term*ne(A[i,perm[i]]))
 det=ne(det+sign*term);terms+=1
print('ROWS',rows,'DET COEFFS',[str(det[i].factor()) if det[i] else '0' for i in range(8)],'terms',terms,'sec',time.monotonic()-st,flush=True)
save({'ring':E,'rows':rows,'determinant':det,'scope':__doc__},str(p/('independent_subcritical_single_RR_three_p_five_minor_'+'_'.join(map(str,rows))+'.sobj')),compress=False)
if det:
 if rows==(2,3,6,21,9):
  m=1+k;q=1+r;u=E.gen()**2/m;expected=E(QQ(16)/42525*k*r*r*q*m**3)*(12*u**3+123*u*u-473*u+468);assert det==expected
  Q=PolynomialRing(QQ,'qq_u');v=Q.gen();F=v**4-2*v**3+2*v*v-v+4;P=12*v**3+123*v*v-473*v+468;res=F.resultant(P);assert res==12087850240;norm=(QQ(16)/42525)**8*res*res*k**8*r**16*q**8*m**24
  save({'ring':E,'rows':rows,'determinant':det,'unit_norm':norm,'QQ_candidate_polynomial':F,'QQ_determinant_polynomial':P,'QQ_resultant':res,'scope':__doc__},str(p/('independent_subcritical_single_RR_three_p_five_minor_'+'_'.join(map(str,rows))+'.sobj')),compress=False)
  print('NORM',norm.factor(),'sec',time.monotonic()-st,flush=True)
 else:
  F=E.modulus();norm=nk(F.resultant(det.lift()));print('NORM',norm.factor(),'sec',time.monotonic()-st,flush=True)
