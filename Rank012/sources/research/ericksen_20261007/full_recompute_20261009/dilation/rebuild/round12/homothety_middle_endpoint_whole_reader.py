"""Units and source data for arbitrary-rate endpoint exclusion."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).absolute().parent;checks={};details={}
for tag,pzero in [('q_zero',False),('p_zero',True)]:
 stem='homothety_middle_endpoint_uniform_quartic'+('_pzero' if pzero else '')
 d=load(str(OUT/(stem+'.sobj')));K=d['field'];S=d['ring'];a,r,A2,B2,kfree,v=K.gens()
 M=d['matrix'];rhs=d['rhs'];rows=d['summary']['rows'];AA=M.matrix_from_rows(rows)
 checks[tag+'_full_original_source_orders']=all(d['checks'].values())
 checks[tag+'_whole_dimensions']=M.dimensions()==(76,4) and len(d['quartic_forcing_free'])==76
 eexpected=-96/r if pzero else -96*r**4
 checks[tag+'_protected_endpoint_minor']=AA.det().subs({v:0})==eexpected
 forcing=d['forcing'];fr=d['forcing_rows'];fc=d['forcing_columns'];FP=forcing.matrix_from_rows_and_columns(fr,fc)
 fexpected=(r*r+a)**2/(a*a*r*r) if pzero else (r*r+a)**2/r**4
 checks[tag+'_forcing_endpoint_minor']=FP.det().subs({v:0})==fexpected
 for i in range(76):
  for j in range(1,4):
   for ob in (A2,B2,kfree):checks[f'{tag}_observation_free_column_{i}_{j}_{ob}']=M[i,j].derivative(ob)==0
  for ob in (A2,B2,kfree):checks[f'{tag}_observation_free_constant_{i}_{ob}']=rhs[i].derivative(ob)==0
 idx=22 if pzero else 21
 coeff=vector(K,M.row(idx))*AA.inverse()
 witness=d['quartic_forcing_free'][idx]-sum(coeff[j]*d['quartic_forcing_free'][rows[j]] for j in range(4))
 checks[tag+'_whole_quartic_source_witness']=witness==d['residuals'][idx]
 positive=20*(1+r*r) if pzero else 20*r*r*(1+r*r)
 checks[tag+'_positive_endpoint_obstruction']=K(witness).subs({v:0})==positive
 for ob in (A2,B2,kfree):checks[tag+'_witness_free_'+str(ob)]=K(witness).derivative(ob)==0
 old=load(str(OUT/('homothety_middle_null_endpoint_cubic'+('_pzero' if pzero else '')+'.sobj')))
 def oldcf(c):return K(c.numerator()(a,r,A2,B2,kfree))/K(c.denominator()(a,r,A2,B2,kfree))
 oldM=old['matrix'].apply_map(oldcf)
 checks[tag+'_cubic_quartic_same_linear_symbol']=M.apply_map(lambda c:c.subs({v:0}))==oldM
 for i,pol in enumerate(d['state_jet']):
  for n,co in enumerate(pol.list()):
   zz=S(co.subs({S.gen(0):0}))
   for ob in (A2,B2,kfree):
    checks[f'{tag}_zero_epsilon_state_observation_independent_{i}_{n}_{ob}']=all(K(c).derivative(ob)==0 for c in zz.coefficients())
 assert all(checks.values()),[kk for kk,vv in checks.items() if not vv]
 details[tag]={'endpoint_minor':str(eexpected),'forcing_endpoint_minor':str(fexpected),
  'endpoint_positive_obstruction':str(positive),'raw_source_checks':sum(d['checks'].values()),
  'plane_parameter':'p=v,q=1' if pzero else 'p=1,q=v'}
out={'status':'PASS_WHOLE_ENDPOINT_RATE_INPUTS','count':len(checks),'details':details,'checks':checks,
 'scope':'Combined with the two uniform analytic estimates in homothety_middle_endpoint_whole_proof.txt, excludes every endpoint approach in epsilon=o(|s|). No assumption on v/s.'}
(OUT/'homothety_middle_endpoint_whole_reader.json').write_text(json.dumps(out,indent=2))
print(json.dumps({kk:vv for kk,vv in out.items() if kk!='checks'}))
