"""All eight original support gauges are relative-regular at the cyclic base.

Only selected p monomials are localized. Every material value remains
symbolic, and all coefficient denominator factors at cyclic are old units.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();sg=load(str(p/'independent_ratio_p_boundary_support_gauges.sobj'));g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));K=g['coefficient_field'];R=g['first_gauged_fields'].base_ring();B=PolynomialRing(QQ,names=('cy_a','cy_r'));a,r=B.gens();F=B.fraction_field();q=1+r;m=1+a*r;ph=K.ring().hom([F(a),a*q/m,1/r,-q/r],F);checks=[];denoms=set();records=[]
def ck(label,v):assert v,label;checks.append(label)
polys=list(g['first_gauged_fields'].list())+list(g['base_fields'].list())
selected_coefficients=[]
for it in sg['items']:
 pol=it['denominator'];ck('selected_p_monomial_%s'%it['chart'],len(pol.dict())==1 and all(sum(map(int,ex[:6]))==3 and not any(ex[6:]) for ex in pol.dict()))
 selected_coefficients+=list(pol.coefficients())
 polys+=list(it['full_field_numerators'].list())+list(it['full_base_numerators'].list())+[it['denominator']]
for pol in polys:
 for co in pol.coefficients():denoms.add(K(co).denominator())
allowed=[F(a),F(r),F(q),F(m),F(a-1)]
for j,co in enumerate(selected_coefficients):
 co=K(co);val=ph(co.numerator())/ph(co.denominator());ck('selected_p_scalar_nonzero_%s'%j,val!=0)
 for ff,k in val.factor():ck('selected_p_scalar_old_unit_%s_%s'%(j,k),any(ff/u in QQ and ff/u for u in allowed))
for j,de in enumerate(denoms):
 val=ph(de);ck('cyclic_den_nonzero_%s'%j,val!=0)
 for ff,k in val.factor():ck('all_material_old_unit_%s_%s'%(j,k),any(ff/u in QQ and ff/u for u in allowed))
 records.append({'original_denominator':de,'cyclic_value':val,'factorization':val.factor()})
save({'source_files':['independent_ratio_p_boundary_support_gauges.sobj','independent_ratio_p_boundary_second_gauge.sobj'],'material_field':F,'cyclic_parameter_map':[a,a*q/m,1/r,-q/r],'denominator_records':records,'checks':checks,'scope':__doc__},str(p/'independent_relative_Newton_cyclic_regularity_reader.sobj'),compress=False)
(p/'independent_relative_Newton_cyclic_regularity_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'distinct_denominators':len(records),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('PASS',len(checks),'DENOMINATORS',len(records),'seconds',time.monotonic()-st,flush=True)
