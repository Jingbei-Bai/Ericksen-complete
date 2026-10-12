"""Independent bridges from unprojected/projected quadratic sources to U5."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).absolute().parent
d=load(str(OUT/'homothety_middle_resonant_whole_coefficients.sobj'))
cb=load(str(OUT/'homothety_middle_null_cubic_face.sobj'))
qd=load(str(OUT/'homothety_axis_plane_raw_quadratic_schur.sobj'))
pd=load(str(OUT/'homothety_axis_plane_projected_quadratics.sobj'))
K=d['field'];a,p,q,r,A2,B2,k1=K.gens();b=a*(1+r*r)/(a+r*r);eta=q+p*r
checks={};denominators={}
def cf(c):return K(c.numerator()(a,b,p,q))/K(c.denominator()(a,b,p,q))
for label,source in [('raw',qd['quadratic_matrix']),('projected',pd['projected_matrix'])]:
 Q2=source.apply_map(cf);M=d[label+'_matrix'];c=d[label+'_constant'];cub=cb['raw_cubic' if label=='raw' else 'forcing_free_cubic']
 for i in range(76):
  checks[f'{label}_null_quadratic_{i}']=Q2[i,0]+r*Q2[i,1]+r*r*Q2[i,4]==0
  checks[f'{label}_d_column_{i}']=M[i,1]==Q2[i,1]+2*r*Q2[i,4]
  checks[f'{label}_h_column_{i}']=M[i,2]==Q2[i,2]+r*Q2[i,5]
  checks[f'{label}_zz_column_{i}']=M[i,3]==Q2[i,9]
  checks[f'{label}_eta_z_column_{i}']=eta*M[i,4]==Q2[i,3]+r*Q2[i,6]
  checks[f'{label}_pure_cubic_{i}']=c[i]==cub[i].constant_coefficient()
  checks[f'{label}_whole_center_{i}']=M[i,4]*d['center_C']+c[i]==0
 allco=list(M.list())+list(c)
 den=K.base_ring().one()
 for z in allco:den=den.lcm(z.denominator())
 dend=K(den).subs({q:-p*r})
 checks[label+'_no_hidden_eta_pole']=dend!=0
 denominators[label]=str(dend.factor())
rows=d['rows'];M=d['projected_matrix'].matrix_from_rows(rows).apply_map(lambda c:c.subs({q:-p*r}))
checks['five_by_five_full_spectrum_unit']=M.det()==4096*r**6*p**7*(r*r+a)/a
checks['nonzero_center']=d['center_C'].subs({q:-p*r})==2*r*r*p**3*(1+r*r)
assert all(checks.values()),[kk for kk,vv in checks.items() if not vv]
out={'status':'PASS_INDEPENDENT_RESONANT_U5_SOURCE_BRIDGES','count':len(checks),'checks':checks,
 'denominator_at_resonance':denominators,'scope':'Checks the whole U5 coefficient source and absence of hidden eta poles; analytic remainder argument is separately audited.'}
(OUT/'homothety_middle_resonant_independent_reader.json').write_text(json.dumps(out,indent=2))
print(json.dumps({kk:vv for kk,vv in out.items() if kk!='checks'}))
