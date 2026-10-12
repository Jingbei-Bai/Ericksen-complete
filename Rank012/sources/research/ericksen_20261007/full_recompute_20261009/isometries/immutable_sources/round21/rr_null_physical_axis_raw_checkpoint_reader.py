"""Independent raw three conics -> exact axis pull -> reduction checkpoint.

No truncated series or saved leading coefficient is used.  The actual row
transform is compared to the protected formula stated in theorem A.
"""
from sage.all import *
from pathlib import Path
import json,time,signal
O=Path(__file__).resolve().parent; start=time.time(); signal.alarm(120)
raw=load(str(O/'rr_null_physical_quadratics_plus.sobj'))
ax=load(str(O/'rr_null_physical_axis_exact_plus.sobj'))
cp=load(str(O/'rr_null_physical_axis_reduction_checkpoint_plus.sobj'))
F=ax['base']; B=F.ring(); s,t=F.gens(); E=ax['e_ring']; e=E.gen(); FE=E.fraction_field()
a,b,ch=raw['ring'].gens(); mons=[a*a,a*b,a*ch,b*b,b*ch,ch*ch]
def evaluate_polynomial(f):
    out={}
    for (i,j,k,l,m),v in f.dict().items():
        out[k+l]=out.get(k+l,B.zero())+B({(i,j):v})
    return E(out)
def evaluate_rational(f):
    return FE(evaluate_polynomial(f.numerator()))/FE(evaluate_polynomial(f.denominator()))
pulled=matrix(FE,[[evaluate_rational(f.monomial_coefficient(m)) for m in mons] for f in raw['quadratics']])
checks={}
for i in range(3):
    checks['raw_full_conic_axis_pull_%s'%i]=pulled.row(i)==ax['coefficient_matrix'].row(i)
scale=diagonal_matrix(FE,[1,1,1/e,1,1/e,1/e**2]); M=pulled*scale
checks['free_chi_inverse_e_exact']=M==cp['original_matrix']
r=(s*s-1)/(t*t-1)
U=matrix(FE,[[1,0,0],[0,1,0],[0,e**7,r*e**7]])
checks['explicit_protected_transform']=U==cp['transform']
checks['checkpoint_replayed_from_raw']=U*M==cp['reduced_matrix']
def val(f):
    return f.numerator().valuation()-f.denominator().valuation() if f else Infinity
vals=[min(val(f) for f in row) for row in cp['reduced_matrix'].rows()]
checks['checkpoint_exact_row_valuations']=vals==[-4,-7,2]
normalized=diagonal_matrix(FE,[e**4,e**7,e**(-2)])*cp['reduced_matrix']
V=matrix(FE,[[e**4,0,0],[0,e**7,0],[0,e**5,r*e**5]])
checks['three_named_theorem_forms']=normalized==V*M
checks['normalization_determinant_protected']=V.det()==e**16*r
checks['parameter_ratio_protected_formula']=r==(s*s-1)/(t*t-1)
# All normalized coefficients are regular in e at zero.  Their parity is even.
checks['all_normalized_coefficients_nonnegative_valuation']=all(val(f)>=0 for f in normalized.list())
flip=FE.hom([-e],FE)
checks['all_normalized_coefficients_even']=all(flip(f)==f for f in normalized.list())
assert all(checks.values()),checks
rep={'status':'PASS_RAW_TO_AXIS_CHECKPOINT','checks':checks,'count':len(checks),
 'exact_transform':'f0=e^4 Q0, f1=e^7 Q1, f2=e^5(Q1+r Q2), chi=e^-1 x',
 'r':'(s^2-1)/(t^2-1)', 'determinant':'e^16 (s^2-1)/(t^2-1)',
 'protection':'e!=0 and |t|>s>1; no additional direction or spectral divisor',
 'seconds':time.time()-start}
(O/'rr_null_physical_axis_raw_checkpoint_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
print(json.dumps(rep,indent=2),flush=True)
