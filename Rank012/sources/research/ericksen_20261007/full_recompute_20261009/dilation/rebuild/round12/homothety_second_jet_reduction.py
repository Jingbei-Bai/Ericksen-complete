"""Original-source tangency of all six eliminated t, with no state pivots."""
from sage.all import *
from pathlib import Path
import time,resource,json,pickle
OUT=Path(__file__).absolute().parent;START=time.time();d=load(str(OUT/'homothety_first_jet_reduction.sobj'))
R=d['ring'];K=R.base_ring();a,b,s,t=K.gens();vs=R.gens();ps=vs[:9];ts=vs[9:15];us=vs[15:];fp=d['scaled_pi_fields'];dn=d['direction_fields'];gam=d['gamma'];sol=d['t_solution'];F=d['fields'];sub=dict(zip(ts,sol))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
def cfD(f,var):return R({mon:cf.derivative(var) for mon,cf in f.dict().items()})
def E(i,f):return sum(f.derivative(p)*fp[i][j] for j,p in enumerate(ps))+cfD(f,s)*dn[i][0]+cfD(f,t)*dn[i][1]
eq=[];labels=[]
for i in range(3):
 for h,z in enumerate(sol):
  eq.append(E(i,z)-F[i][9+h].subs(sub)-2*gam[i]*z);labels.append((i,h))
 print('tangent direction',i,'seconds',round(time.time()-START,2),flush=True)
M=matrix(K,[[z.monomial_coefficient(u) for u in us] for z in eq]);rows=list(M.transpose().pivots());MM=M.matrix_from_rows(rows)
assert len(rows)==2
rhs=vector(R,[-eq[k].subs({u:0 for u in us}) for k in rows]);usol=MM.inverse()*rhs
res=[z.subs(dict(zip(us,usol))) for z in eq]
orig=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
C3=[R({ex:K(cf) for ex,cf in z.items()}).subs(sub) for z in orig['C3']]
print('u rows',rows,'det',MM.det().factor(),'residuals',sum(bool(x) for x in res),'terms',[len(z.dict()) for z in res if z],flush=True)
save({'source':d,'ring':R,'u_solution':usol,'u_matrix':M,'u_pivot_rows':rows,'u_pivot_det':MM.det(),'t_tangency_residuals':res,'t_tangency_labels':labels,'C3_substituted':C3},str(OUT/'homothety_second_jet_reduction.sobj'))
(OUT/'homothety_second_jet_reduction.json').write_text(json.dumps({'status':'EXACT_NECESSARY_SOURCE_TANGENCY','u_rows':rows,'u_det':str(MM.det().factor()),'tangency_nonzero_count':sum(bool(x) for x in res),'terms':[len(z.dict()) for z in res if z],'C3_terms':[len(z.dict()) for z in C3],'seconds':time.time()-START},indent=2))
print('COMPLETE',time.time()-START,flush=True)
