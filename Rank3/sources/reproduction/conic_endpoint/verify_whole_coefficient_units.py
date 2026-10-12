"""Census all generic 28-source coefficients, including positive Newton layers."""
from sage.all import *
from pathlib import Path
import hashlib,json,time
out=Path(__file__).absolute().parent;p=out/'source_inputs'
start=time.monotonic()
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
old=full['coefficient_field'];OK=qf['ring'].base_ring()
K=PolynomialRing(QQ,names=('k','r')).fraction_field();k,r=K.gens();m=1+k;q=1+r
ph=old.ring().hom([k/r,k*q/(r*m),1/r,-q/r],K)
po=OK.ring().hom([k/r,r,1,-1],K)
denoms={'original':set(),'normalized':set()};counts={}
def add(label,coefficients,F,which):
    n=0
    for c in coefficients:denoms[which].add(F(c).denominator());n+=1
    counts[label]=counts.get(label,0)+n
for i,f in enumerate(qf['full_R_sources']):add('R'+str(i),f.coefficients(),OK,'normalized')
add('Q0',qf['full_Q0'].coefficients(),OK,'normalized')
for name in ['D','C','L','Q0_L_cofactors']:
    for f in qf[name].list():add(name,f.coefficients(),OK,'normalized')
for item in ev['items']:
    label='even'+str(item['Q_index'])
    add(label,item['full_source'].coefficients(),OK,'normalized')
    add(label+' square',[item['square_unit']],OK,'normalized')
    for f in item['A_p_quadratic']:add(label+' A',f.coefficients(),OK,'normalized')
refs=set()
for record in full['records']:
    for sign,li,ri in record['full_signed_product_DAG']:refs.update([li,ri])
for ref in refs:add('H full product pool',full['polynomial_pool'][ref].values(),old,'original')
for item in hs['items']:
    for dd in item['L_cofactors_original']:add('H L cofactors',dd.values(),old,'original')
records=[];factors=set()
for which,items in denoms.items():
    hom=ph if which=='original' else po
    for dd in items:
        value=K(hom(dd));assert value,('undefined coefficient at cyclic base',str(dd))
        local_factors=[]
        for poly in [value.numerator(),value.denominator()]:
            for fac,power in poly.factor():
                fac=K(fac)
                assert any((fac/u).numerator().total_degree()==0 and
                           (fac/u).denominator().total_degree()==0
                           for u in [K.one(),k,r,m,q,k-r]),('unprotected factor',str(fac))
                factors.add(str(fac));local_factors.append(str(fac))
        records.append({'field':which,'source_denominator':str(dd),'cyclic_value':str(value),'factors':local_factors})
inputs=[
'independent_delta_subcritical_Q0_full_source.sobj',
'independent_delta_subcritical_H_initials.sobj',
'independent_ratio_global_R_derivative_degreefall.sobj',
'independent_delta_subcritical_even_full_source.sobj']
files=[]
for name in inputs:
    path=p/name;digest=hashlib.file_digest(path.open('rb'),'sha256').hexdigest()
    files.append({'name':name,'bytes':path.stat().st_size,'sha256':digest})
result={'status':'PASS','scope':__doc__,'coefficients_by_group':counts,
 'total_input_coefficient_occurrences':sum(counts.values()),
 'distinct_denominators':{key:len(val) for key,val in denoms.items()},
 'cyclic_material_factors':sorted(factors),'source_denominator_records':records,
 'exact_inputs':files,'total_input_bytes':sum(item['bytes'] for item in files),
 'seconds':time.monotonic()-start}
(out/'whole_coefficient_units.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS all coefficients',sum(counts.values()),'distinct denominators',len(records),'factors',sorted(factors),flush=True)
print('INPUT BYTES',result['total_input_bytes'],'seconds',result['seconds'],flush=True)
