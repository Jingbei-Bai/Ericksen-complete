"""Full-source quadratic cuts the actual constant-Z infinity amplitude line."""
from pathlib import Path
import time,json
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'general_curve_constant_Z_quadratic_relations.sobj'));R=d['ring'];K=R.base_ring();a,r,H,Z=K.gens();C=K.ring();xs=R.gens()
f=d['quadratics'][0];bad=R({m:c for m,c in f.dict().items() if sum(m)==2 and all(m[j]==0 for j in (0,2,3,5))})
c=f.monomial_coefficient(xs[1]**2);unit=Z*(r+1)*(Z-1)**2/2;P=C(c/unit)
checks={'quadratic_factorization':bad==c*(xs[1]**2+(r+1)*Z**2/r*xs[4]**2),'P_independent_a':not P.derivative(C.gen(0)),'P_independent_r':not P.derivative(C.gen(1))}
U=PolynomialRing(QQ,'z');z=U.gen();coeffs=[U.zero()]*(P.degree(C.gen(2))+1)
for mon,co in P.dict().items():coeffs[mon[2]]+=co*z**int(mon[3])
g=coeffs[0];multipliers=[U.one()]+[U.zero()]*(len(coeffs)-1)
for j in range(1,len(coeffs)):
 gg,s,t=g.xgcd(coeffs[j]);multipliers=[v*s for v in multipliers];multipliers[j]+=t;g=gg
checks['Bezout_identity']=sum(m*p for m,p in zip(multipliers,coeffs))==g
checks['coefficient_gcd_has_only_Z_old_guards']=all(z*(z-1)%fac==0 for fac,e in g.factor())
assert all(checks.values())
save({'ring':R,'full_quadratic':f,'bad_line_restriction':bad,'factor_coefficient':c,'protected_unit':unit,'parameter_polynomial':P,'H_coefficients':coeffs,'Bezout_multipliers':multipliers,'coefficient_gcd':g,'checks':checks,'scope':'On every allowed fixed spectrum and every allowed constant Z, a full necessary quadratic restricts nontrivially to the previously surviving two-pair amplitude P1 for generic H. Its projective zero set has at most two points. Other supports and total fiber dimension not concluded.'},str(out/'general_curve_constant_Z_boundary_cut.sobj'))
report={'checks':checks,'coefficient_gcd':str(g.factor()),'coefficient_degrees':[int(p.degree()) for p in coeffs],'Bezout_degrees':[int(p.degree()) if p else -1 for p in multipliers],'elapsed':time.monotonic()-start}
(out/'general_curve_constant_Z_boundary_cut.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
