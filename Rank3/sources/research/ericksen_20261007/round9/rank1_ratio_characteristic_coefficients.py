"""Compact coefficient table for the affine-pair determinant strata.

Use r=q-1, c=q*T/r, h=H/r, b=a*q/(1+a*r).
T excludes 0,1; H, H+T, H+1, 1-H-T exclude zero.
"""
from pathlib import Path
import time,json
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'rank1_ratio_characteristic_linear_pair.sobj'));O=d['ring'];u,v,x,s,y=O.gens()
f,g=d['linear_source_rows'];j=d['remaining_mixed_source']
old_coeff={'lam':f.monomial_coefficient(v*v*y),'A':f.monomial_coefficient(v*s*s),
 'B':f.monomial_coefficient(y*s*s),'C':f.monomial_coefficient(u*s),'D':f.monomial_coefficient(x*s),
 'E':g.monomial_coefficient(u*s*s),'F':g.monomial_coefficient(x*s*s),'G':g.monomial_coefficient(v*s),
 'Hcoef':g.monomial_coefficient(y*s),'C0':f.monomial_coefficient(v),'Jxy':j.monomial_coefficient(x*y*s)}
CR=PolynomialRing(QQ,names=('aa','qq','TT','HH'));K=CR.fraction_field();a,q,T,H=map(K,CR.gens());r=q-1
hom=O.base_ring().ring().hom([a,a*q/(1+a*r),q*T/r,H/r],K)
coeff={k:hom(co.numerator())/hom(co.denominator()) for k,co in old_coeff.items()}
coeff['alpha']=coeff['Hcoef']+coeff['lam']*coeff['D']
coeff['E_minus_A']=coeff['E']-coeff['A']
records={k:str(c.factor()) for k,c in coeff.items()}
save({'coefficient_ring':CR,'field':K,'parameters':(a,q,T,H),'coefficients':coeff,'old_coefficients':old_coeff,
 'records':records,'source_file':'rank1_ratio_characteristic_linear_pair.sobj'},str(out/'rank1_ratio_characteristic_coefficients.sobj'))
(out/'rank1_ratio_characteristic_coefficients.json').write_text(json.dumps({'records':records,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
print(json.dumps(records,indent=2),flush=True)
