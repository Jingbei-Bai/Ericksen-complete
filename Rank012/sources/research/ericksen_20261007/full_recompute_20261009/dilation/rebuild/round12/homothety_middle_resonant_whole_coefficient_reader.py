"""Exact coefficients for the whole resonant estimate; no rate sampling."""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).absolute().parent;START=time.time()
d=load(str(OUT/'homothety_middle_null_cubic_face.sobj'));qd=load(str(OUT/'homothety_axis_plane_raw_quadratic_schur.sobj'))
pd=load(str(OUT/'homothety_axis_plane_projected_quadratics.sobj'))
K=d['parameter_field'];a,p,q,r,A2,B2,k1=K.gens();b=a*(1+r*r)/(a+r*r);eta=p*r+q
S=d['ring'];E,t2,h2,z2=S.gens();C=-2*r*p*p*q*(1+r*r)
def cf4(c):return K(c.numerator()(a,b,p,q))/K(c.denominator()(a,b,p,q))
checks={};out={}
for name,cub,Q2 in [('projected',d['forcing_free_cubic'],pd['projected_matrix']),('raw',d['raw_cubic'],qd['quadratic_matrix'])]:
 MM=[];const=[]
 for i,f in enumerate(cub):
  vals=[K(f.monomial_coefficient(v)) for v in (E,t2,h2,z2)]
  checks[f'{name}_affine_cubic_{i}']=f==sum(vals[j]*S.gen(j) for j in range(4))+f.constant_coefficient()
  checks[f'{name}_resonant_z_column_{i}']=vals[3].subs({q:-p*r})==0
  row=vals[:3]+[cf4(Q2[i,9]),vals[3]/eta];MM.append(row);const.append(K(f.constant_coefficient()))
 mat=matrix(K,MM);cv=vector(K,const)
 checks[f'{name}_moving_cubic_center']=mat*vector(K,[0,0,0,0,C])+cv==0
 for j in (1,2,3,4):
  for i in range(mat.nrows()):
   checks[f'{name}_observation_free_column_{i}_{j}']=all(mat[i,j].derivative(v)==0 for v in (A2,B2,k1))
 checks[f'{name}_constant_observation_free']=all(z.derivative(v)==0 for z in cv for v in (A2,B2,k1))
 out[name+'_matrix']=mat;out[name+'_constant']=cv
rows=[46,20,12,21,34];limit=out['projected_matrix'].apply_map(lambda c:c.subs({q:-p*r}));block=limit.matrix_from_rows(rows)
expected=4096*r**6*p**7*(r*r+a)/a
checks['resonant_five_by_five_unit_determinant']=block.det()==expected
checks['nonzero_resonant_center']=C.subs({q:-p*r})==2*r*r*p**3*(1+r*r)
# The sign of delta is not constrained: X=z^2/delta is allowed to be signed.
# At delta=-tau^2 the four E,t2,h2,Y columns and constant change sign, while
# the X column does not. Verify the same algebraic contradiction on that side.
sg=diagonal_matrix(K,[-1,-1,-1,1,-1]);negative=block*sg
cc=vector(K,[out['projected_constant'][i].subs({q:-p*r}) for i in rows])
negative_solution=negative.inverse()*cc
checks['negative_delta_ramified_solution']=negative_solution==vector(K,[0,0,0,0,2*r*r*p**3*(1+r*r)])
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'field':K,'eta':eta,'center_C':C,'rows':rows,'limiting_matrix':block,
      'determinant':expected,'checks':checks,**out},str(OUT/'homothety_middle_resonant_whole_coefficients.sobj'))
record={'status':'PASS','checks':checks,'count':len(checks),'determinant':str(expected.factor()),
 'center_C':str(C.factor()),'center_at_resonance':str(C.subs({q:-p*r}).factor()),
 'scope':'Whole analytic Taylor coefficient identities; proof of the rate estimate is in the companion text. Both signed delta sides are retained.',
 'seconds':float(time.time()-START)}
(OUT/'homothety_middle_resonant_whole_coefficient_reader.json').write_text(json.dumps(record,indent=2))
print(json.dumps({'status':'PASS','count':len(checks),'seconds':record['seconds']}),flush=True)
