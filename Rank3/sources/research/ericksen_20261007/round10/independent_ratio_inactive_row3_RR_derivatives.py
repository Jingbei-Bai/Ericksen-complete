"""All-parameter ordinary derivative source on p31=p32=0, no p localization."""
from pathlib import Path
from sage.all import load,save,PolynomialRing
import time,json
out=Path(__file__).resolve().parent;start=time.monotonic()
raw=load(str(out.parent/'round8/rank3_joint_fiber_leading.sobj'))
gauge=load(str(out/'independent_ratio_p_boundary_second_gauge.sobj'))
g=load(str(out/'independent_ratio_generic_canonical_normal_28.sobj'))
K=g['coefficient_field'];a,r,H,Z=K.gens();oldK=raw['ring'].base_ring();ph=oldK.ring().hom(g['parameter_map'],K);cache={}
R=PolynomialRing(K,names=tuple('ir_p%d'%i for i in range(6))+tuple('ir_rr%d'%i for i in range(3)));v=R.gens()
def cv(c):
 c=oldK(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def cp(f,degree=None):return R({tuple(map(int,e)):cv(c) for e,c in f.dict().items() if int(e[4])==0 and int(e[5])==0 and (degree is None or sum(map(int,e))==degree)})
items=[];checks=0
field=cp(gauge['first_gauged_fields'][2,8],3)
unit=K(r*(H+1)/(H*(1+r)));Z0=field/unit
assert field.derivative(v[8])==0;checks+=1
print('FIELD3',len(field.dict()),'Z0 terms',len(Z0.dict()),flush=True)
for idx in (0,1,5,6):
 f=raw['quartic_source_rows'][idx];degree=max(sum(map(int,e)) for e in f.dict());top=cp(f,degree)
 RRdeg=int(top.degree(v[8]));ders=[]
 for k in range(1,RRdeg+1):
  der=top.derivative(v[8],k);factor=der.factor() if der else None
  ders.append({'order':k,'derivative':der,'factorization':factor,'principal_iterate':field**k*der})
  print('Q',idx,'dRR3',k,'terms',len(der.dict()),'factor',factor,flush=True)
 checks+=1;items.append({'source_index':idx,'original_ordinary_degree':int(degree),'restricted_initial':top,'RR3_degree':RRdeg,'derivatives':ders})
save({'ring':R,'coefficient_field':K,'parameters':K.gens(),'original_parameter_map':g['parameter_map'],'row_constraint':(v[4],v[5]),'field_ordinary_cubic':field,'source_unit':unit,'missing_Z0':Z0,'items':items,'checks':checks,'scope':__doc__,'logical_qualification':'Ordinary degree+2 principal field is field_ordinary_cubic*d/dRR3. Its iterated leading components are necessary only through differentiation of the genuine full original sources; restriction to the inactive row is evaluation at initial points, not assumed flow-invariance of that row.'},str(out/'independent_ratio_inactive_row3_RR_derivatives.sobj'),compress=False)
(out/'independent_ratio_inactive_row3_RR_derivatives.json').write_text(json.dumps({'checks':checks,'field_terms':len(field.dict()),'sources':[{'source_index':z['source_index'],'RR3_degree':z['RR3_degree'],'derivative_terms':[len(t['derivative'].dict()) for t in z['derivatives']]} for z in items],'seconds':time.monotonic()-start},indent=2)+'\n');print('DONE',checks,time.monotonic()-start,flush=True)
