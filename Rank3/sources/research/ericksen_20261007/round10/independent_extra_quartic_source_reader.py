"""Independent provenance audit for six additional necessary cubic sources."""
from pathlib import Path
from functools import lru_cache
import json,time
from sage.all import load,save,prod

out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic();checks={}
d=load(str(out/'general_curve_extra_quartic_cubics.sobj'))
normal=load(str(out/'rank3_quartic_diagonal_velocity_form.sobj'))
joint=load(str(out/'joint_fields_direct.sobj'))
quartics=load(str(base/'round8/rank3_joint_fiber_leading.sobj'))
curve=load(str(out/'rank1_general_slope_quartics.sobj'))
full=load(str(out/'general_curve_dimensionless_full.sobj'))
OR=joint['ring'];ox=OR.gens()
for pos,(qi,i,xi,sgn) in enumerate(normal['row_specs']):
    Wi=joint['fields'][i,6+i];k=normal['protected_coefficients'][i]
    checks['original_quartic_'+str(qi)]=OR(quartics['quartic_source_rows'][qi])==normal['remainders'][pos]+sgn*ox[xi]*Wi/k
    checks['original_diagonal_field_'+str(i)]=normal['diagonal_velocities'][i]==Wi
print('original six quartics replayed',round(time.monotonic()-start,2),flush=True)

R=d['ring'];K=R.base_ring();a,r,H,Z,S,N=K.gens();q=1+r;b=a*q/(1+a*r)
xs=R.gens();w1,m1,w2,m2,w3,m3=xs
pi=[w3,w2,((1+H)*S-Z)*m3/(2*Z)+H*w3,w1,-(1+H*S-Z)*m2/(2*H)-Z*w2,-(1-S)*m1/(2*H)-Z*w1/H]
rr=[(Z-H*S)*m1/(H*Z),S*m2/Z,-m3/H]
images=pi+rr
checks['full_coordinate_images']=images==d['images']
@lru_cache(None)
def cc(c):return K(c.numerator()(a,b,H/r,q*Z/r))/K(c.denominator()(a,b,H/r,q*Z/r))
def cc5(c):return K(c.numerator()(a,b,H/r,q*Z/r,q*S))/K(c.denominator()(a,b,H/r,q*Z/r,q*S))
def ev(poly,im=images,coeff=cc):
    return sum((coeff(c)*prod([im[j]**int(e) for j,e in enumerate(mm) if e]) for mm,c in poly.dict().items()),R.zero())
allowed=[a,r,q,1+a*r,H,Z,H-Z,1-Z,1+H,1+H-Z,1-S,1+H*S-Z,(1+H)*S-Z]
def is_unit(c):
    c=K(c)
    if not c:return False
    for p in [c.numerator(),c.denominator()]:
        for f,e in p.factor():
            if f.is_constant():continue
            if not any(all((K(f)/g).derivative(x)==0 for x in K.gens()) for g in allowed):return False
    return True
Wi=[];Ti=[];piv=[];ks=[]
for i in range(3):
    checks['ratio_graph_'+str(i)]=ev(curve['ratio_graph'][i],pi,cc5)==rr[i]
    piv.append(cc5(curve['pivots'][i]));ks.append(cc(normal['protected_coefficients'][i]))
    checks['pivots_protected_'+str(i)]=is_unit(piv[i]) and is_unit(ks[i])
    wi=ev(joint['fields'][i,6+i]);Wi.append(wi)
    ti=full['rows'][4+i]-piv[i]*wi;Ti.append(ti)
    checks['actual_diagonal_source_'+str(i)]=wi==d['diagonal_velocities'][i] and ti==d['diagonal_lower_remainders'][i]
    checks['lower_degree_'+str(i)]=ti.total_degree()<=2

for pos,(qi,i,xi,sgn) in enumerate(normal['row_specs']):
    rem=ev(normal['remainders'][pos])
    f=rem-sgn*pi[xi]*Ti[i]/(ks[i]*piv[i])
    checks['new_source_'+str(qi)]=f==d['rows'][pos]
    checks['linear_source_identity_'+str(qi)]=f+sgn*pi[xi]*full['rows'][4+i]/(ks[i]*piv[i])==rem+sgn*pi[xi]*Wi[i]/ks[i]
    checks['new_degree_'+str(qi)]=f.total_degree()<=3
    top=R({tuple(mm):c for mm,c in f.dict().items() if sum(mm)==3})
    checks['new_highest_'+str(qi)]=top==d['leading_rows'][pos]
    print('new Q',qi,'replayed',round(time.monotonic()-start,2),flush=True)
checks['all_six_rows']=len(d['rows'])==6 and len(d['leading_rows'])==6 and d['row_specs']==normal['row_specs']
assert all(checks.values()),[key for key,val in checks.items() if not val]
payload={'checks':checks,'row_specs':normal['row_specs'],'source_files':[
 'round8/rank3_joint_fiber_leading.sobj','joint_fields_direct.sobj','rank3_quartic_diagonal_velocity_form.sobj',
 'rank1_general_slope_quartics.sobj','general_curve_dimensionless_full.sobj','general_curve_extra_quartic_cubics.sobj'],
 'scope':'Six genuine cubic necessary consequences of full original quartics and actual diagonal curve tangencies. All a,r retained and only old noncharacteristic/main guards divided. No endpoint or dimension conclusion follows from provenance alone.'}
save(payload,str(out/'independent_extra_quartic_source_reader.sobj'))
report={'status':'EXACT_ORIGINAL_SIX_QUARTIC_SOURCE_PASS','checks':checks,'check_count':len(checks),'elapsed_seconds':time.monotonic()-start}
(out/'independent_extra_quartic_source_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({key:val for key,val in report.items() if key!='checks'}),flush=True)
