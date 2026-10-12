"""Linear point-state chart retaining both complete angular vectors and pressures."""
from sage.all import *
from pathlib import Path
import time,json
OUT=Path(__file__).absolute().parent;START=time.time();d=load(str(OUT/'homothety_first_jet_reduction.sobj'));R=d['ring'];K=R.base_ring();a,b,s,t=K.gens();P=R.gens()[:9];pi=d['pi'];n=vector(R,d['direction']);L=diagonal_matrix(K,[1,1/a,1/b]);Bvec=(2*pi.transpose()*L-(L*pi).trace()*matrix.identity(R,3))*n/2
obs=list(d['source_rotation_vector'])+list(Bvec)+list(d['longitudinal_pressures'][:2]);M=matrix(K,[[z.monomial_coefficient(p) for p in P] for z in obs]);rk=M.rank();cols=list(M.pivots());print('rank',rk,'cols',cols,flush=True)
N=matrix(K.ring(),8,9,[K.ring()(2*a*b*z) for z in M.list()]);cof=[]
if rk==8:
 for free in range(9):
  idx=[z for z in range(9) if z!=free];dd=N.matrix_from_columns(idx).det()/K(2*a*b)**8
  if dd:print('free',free,'det',dd.factor(),flush=True)
  cof.append((-1)**free*K(dd))
ker=vector(K,[z/cof[2] for z in cof]);assert M*ker==vector(K,[0]*8);assert ker[2]==1
lead=[]
for z in d['residuals']:
 if not z:continue
 zz=sum(co*prod(ker[i]**m[i] for i in range(9)) for m,co in z.dict().items() if sum(m)==2)
 lead.append(K(zz));print('quadratic_fiber_leader',K(zz).factor(),flush=True)
save({'source':d,'matrix':M,'observables':obs,'rank':rk,'target_angular_scaled_vector':Bvec,'fiber_kernel':ker,'quadratic_fiber_leaders':lead},str(OUT/'homothety_invariant_linear_chart.sobj'))
print('COMPLETE',round(time.time()-START,2),flush=True)
