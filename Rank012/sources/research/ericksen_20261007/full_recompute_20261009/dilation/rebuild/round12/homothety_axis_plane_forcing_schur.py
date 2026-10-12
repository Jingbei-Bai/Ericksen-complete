"""All original first normal equations, arbitrary angular/pressure forcing."""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).absolute().parent;START=time.time()
d=load(str(OUT/'homothety_axis_plane_normal_sources.sobj'));K=d['parameter_field'];a,b,p,q=K.gens()
R=PolynomialRing(K,names=('A0','A1','A2','B0','B1','B2','k1','k2','h','z','s','t','e'))
ob=R.gens()[:8];h,z,s,t,e=R.gens()[8:]
N=d['state_matrix'];O=d['observation_matrix'];M=O.stack(N);rows=d['combined_pivot_rows'];cols=d['combined_pivot_columns']
free={2:R(0),5:R(0),8:h,9:z};rhs=[]
# Source and weighted-target angular maps at pi=p E13+q E23.
oforce=[0,0,p*s+q*t,0,0,p*s+q*t/a,q*s+p*t,-q*s-p*t]
rhs.extend([R(oforce[i])-e*ob[i] for i in range(8)])
forcing=d['matrix'][:,17:]
rhs.extend([R(forcing[i,0])*e+R(forcing[i,1])*s+R(forcing[i,2])*t for i in range(N.nrows())])
A=M.matrix_from_rows_and_columns(rows,cols)
v=vector(R,[-rhs[i]-sum(M[i,j]*val for j,val in free.items()) for i in rows]);sol=A.inverse()*v
state=[R(0)]*17
for j,val in free.items():state[j]=val
for j,val in zip(cols,sol):state[j]=val
res=list(M*vector(R,state)+vector(R,rhs))
assert all(res[i]==0 for i in rows)
for rr in res:
 assert rr.subs({e:0})==0
H=[R(rr/e) for rr in res]
CM=matrix(K,[[rr.monomial_coefficient(w) for w in ob]+[rr.constant_coefficient()] for rr in H])
assert all(rr==sum(CM[i,j]*ob[j] for j in range(8))+CM[i,8] for i,rr in enumerate(H))
cr=list(CM.transpose().pivots());cc=list(CM.pivots())
minor=CM.matrix_from_rows_and_columns(cr,cc).det()
strata={}
for name,subst in [('p_nonzero_q_zero',{q:0}),('p_zero_q_nonzero',{p:0})]:
 mm=CM.apply_map(lambda v:v.subs(subst));rr=list(mm.transpose().pivots());cl=list(mm.pivots())
 strata[name]={'rank':len(cl),'rows':rr,'cols':cl,'minor':str(mm.matrix_from_rows_and_columns(rr,cl).det().factor())}
summary={'status':'FIRST_NORMAL_FORCING_ONLY_NOT_INCOMING_ARC_EXCLUSION','rank':len(cc),
 'rows':cr,'columns':cc,'minor':str(minor.factor()),'strata':strata,
 'residuals':[(i,str(rr.factor()) if rr else '0') for i,rr in enumerate(H) if rr],
 'kernel':[[str(z.factor()) if z else '0' for z in row] for row in CM.right_kernel().basis_matrix().rows()],
 'seconds':float(time.time()-START)}
save({'ring':R,'state_linear_graph':state,'linear_residuals':res,'forcing_matrix':CM,
      'source':d,'summary':summary},str(OUT/'homothety_axis_plane_forcing_schur.sobj'))
(OUT/'homothety_axis_plane_forcing_schur.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary),flush=True)
