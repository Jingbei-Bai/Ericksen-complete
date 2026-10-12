"""Only the axis source used in the proof, from freshly regenerated full conics.

The historical exploratory asymptotic script continues to an unused 6x6
series determinant.  This file performs the exact 18 coefficient pull and
the protected row operation actually used in the current theorem.
"""
from sage.all import *
from pathlib import Path
import json,time,hashlib
O=Path(__file__).absolute().parent;src=O/'replay/round21';st=time.time()
raw=load(str(src/'rr_null_physical_quadratics_plus.sobj'))
K=raw['field'];QI=K.ring().base_ring()
B=PolynomialRing(QI,names=('s','t'));F=B.fraction_field();s,t=F.gens()
E=PolynomialRing(F,'e');e=E.gen();FE=E.fraction_field()
def pol(z):
    dd={}
    for (i,j,k,l,m),co in z.dict().items():dd[k+l]=dd.get(k+l,B(0))+B({(i,j):co})
    return E(dd)
def ev(z):return FE(pol(z.numerator()))/FE(pol(z.denominator()))
a,b,ch=raw['ring'].gens();mm=[a*a,a*b,a*ch,b*b,b*ch,ch*ch]
M0=matrix(FE,[[ev(f.monomial_coefficient(m)) for m in mm] for f in raw['quadratics']])
M=M0*diagonal_matrix(FE,[1,1,1/e,1,1/e,1/e**2])
r=(s*s-1)/(t*t-1)
U=matrix(FE,[[1,0,0],[0,1,0],[0,e**7,r*e**7]])
N=U*M
def val(f):return f.numerator().valuation()-f.denominator().valuation() if f else Infinity
checks={'row_valuations':[min(val(x) for x in row) for row in N.rows()]==[-4,-7,2]}
V=matrix(FE,[[e**4,0,0],[0,e**7,0],[0,e**5,r*e**5]])
normalized=V*M
checks['explicit_transform']=normalized==diagonal_matrix(FE,[e**4,e**7,e**(-2)])*N
checks['protected_transform_determinant']=V.det()==e**16*r
checks['all_coefficients_regular']=all(val(x)>=0 for x in normalized.list())
checks['all_coefficients_even']=all(x(-e)==x for x in normalized.list())
exact_path=src/'rr_null_physical_axis_exact_plus.sobj'
checkpoint_path=src/'rr_null_physical_axis_reduction_checkpoint_plus.sobj'
if exact_path.exists():
 old=load(str(exact_path));checks['all18_axis_coefficients']=M0==old['coefficient_matrix']
if checkpoint_path.exists():
 cp=load(str(checkpoint_path));checks['complete_checkpoint']=N==cp['reduced_matrix'] and U==cp['transform']
assert all(checks.values()),checks
save({'base':F,'e_ring':E,'coefficient_matrix':M0,'tag':'plus'},str(exact_path))
save({'field':F,'e_ring':E,'original_matrix':M,'reduced_matrix':N,'transform':U,'history':[]},str(checkpoint_path))
save({'field':F,'e_ring':E,'original_matrix':M,'reduced_matrix':N,'transform':U,'normalized':normalized,'checks':checks},str(O/'axis_pull_fresh.sobj'))
report={'status':'PASS','checks':checks,'count':len(checks),'seconds':time.time()-st,'fresh_source_sha256':hashlib.sha256((src/'rr_null_physical_quadratics_plus.sobj').read_bytes()).hexdigest(),'unused_computation':'No six-conic Jacobian determinant series is required by the present theorem.'}
(O/'axis_pull_fresh.json').write_text(json.dumps(report,indent=2),encoding='utf8');print(json.dumps(report,indent=2),flush=True)
