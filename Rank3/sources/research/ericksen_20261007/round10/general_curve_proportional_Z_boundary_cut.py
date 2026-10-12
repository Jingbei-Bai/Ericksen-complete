"""Uniform full-source cut of the Z=cH two-pair infinity amplitude line."""
from pathlib import Path
import time,json
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'general_curve_proportional_Z_quadratic_relations.sobj'));R=d['ring'];K=R.base_ring();a,r,H,c=K.gens();C=K.ring();xs=R.gens()
f=d['quadratics'][0];bad=R({m:co for m,co in f.dict().items() if sum(m)==2 and all(m[j]==0 for j in (0,1,2,4))})
coef=f.monomial_coefficient(xs[3]**2);unit=4*c*(r+1)*(c-1)**2*H**2/r;P=C(coef/unit)
checks={'quadratic_factorization':bad==coef*(xs[3]**2+(r+1)*c**2*xs[5]**2),
 'P_independent_a':not P.derivative(C.gen(0)),'P_independent_r':not P.derivative(C.gen(1)),
 'H_constant_coefficient':P.subs({C.gen(2):0})==32*(1-C.gen(3))}
assert all(checks.values())
save({'ring':R,'full_quadratic':f,'bad_line_restriction':bad,'factor_coefficient':coef,'protected_unit':unit,'parameter_polynomial':P,'checks':checks,'scope':'For every allowed spectrum and every c !=0,1, a full necessary quadratic cuts the previously surviving Z=cH two-pair amplitude P1 to at most two points for generic H. This includes c=-1 despite the intermediate kernel pivot containing c+1: the final source multipliers and identity are denominator-cleared.'},str(out/'general_curve_proportional_Z_boundary_cut.sobj'))
report={'checks':checks,'H_constant_coefficient':'32*(1-c)','elapsed':time.monotonic()-start}
(out/'general_curve_proportional_Z_boundary_cut.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
