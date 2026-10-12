"""Full-spectrum leading source at an angular-observation infinity axis.

This is a necessary boundary calculation, not an exclusion of incoming arcs.
The scale-one term has lower pi degree and is intentionally omitted only
in this explicitly labelled leading system.
"""
from sage.all import *
from pathlib import Path
import pickle,json,time

OUT=Path(__file__).absolute().parent
START=time.time()
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
K=FractionField(PolynomialRing(QQ,names=('a','b')))
a,b=K.gens()
R=PolynomialRing(K,names=('h','p','q')+tuple('t'+str(i) for i in range(6))+('v0','v1'))
h,p,q=R.gens()[:3];ts=R.gens()[3:9];us=R.gens()[9:]
ps=[h,0,p,0,a*(b-1)/(a-b)*h,q,0,0,b*(a-1)/(a-b)*h]
vals=list(map(R,ps))+list(ts)+list(us)
def pull(comp):
    return sum(R(K(cf))*prod(vals[i]**ex[i] for i in range(17)) for ex,cf in comp.items())
eq=vector(R,[pull(comp) for comp in src['Fraw'][2][:9]])
M=matrix(K,[[z.monomial_coefficient(tt) for tt in ts] for z in eq])
rows=list(M.transpose().pivots());cols=list(M.pivots())
remaining=[j for j in range(6) if j not in cols]
rhs=vector(R,[-eq[i].subs({ts[j]:0 for j in cols}) for i in rows])
sol=M.matrix_from_rows_and_columns(rows,cols).inverse()*rhs
sub={ts[j]:sol[i] for i,j in enumerate(cols)}
res=[z.subs(sub) for z in eq]
summary={'status':'NECESSARY_AXIS_LEADING_ONLY','symbol_rank':M.rank(),
 'rows':rows,'columns':cols,'free_t_indices':remaining,
 't_graph':{str(tt):str(v.factor()) if v else '0' for tt,v in sub.items()},
 'first_residuals':[str(z.factor()) for z in res if z],
 'source_keys':list(src)}
print(json.dumps(summary,indent=2),flush=True)
plane_sub={h:0,**sub}
plane_c3=[pull(comp).subs(plane_sub) for comp in src['C3']]
plane_f3t=[pull(comp).subs(plane_sub) for comp in src['Fraw'][2][9:15]]
plane_q4=[pull(comp).subs(plane_sub) for comp in src['Q4']]
print('PLANE C3',[str(z.factor()) for z in plane_c3 if z],flush=True)
print('PLANE F3t',[str(z.factor()) for z in plane_f3t if z],flush=True)
print('PLANE Q4',[str(z.factor()) for z in plane_q4 if z],flush=True)
summary['plane_C3']=[str(z.factor()) for z in plane_c3 if z]
summary['plane_F3t']=[str(z.factor()) for z in plane_f3t if z]
summary['plane_Q4']=[str(z.factor()) for z in plane_q4 if z]
exact_plane={**plane_sub,ts[0]:p*q,us[0]:0,us[1]:0}
checks={}
for i,comp in enumerate(src['Fraw'][2]):
    checks['plane_F3_'+str(i)]=not pull(comp).subs(exact_plane)
for family in ['C3','Q4','basis']:
    seq=src[family][:27] if family=='basis' else src[family]
    for i,comp in enumerate(seq):
        checks['plane_'+family+'_'+str(i)]=not pull(comp).subs(exact_plane)
assert all(checks.values())
summary['exact_boundary_plane_checks']=len(checks)
summary['exact_boundary_plane_all_pass']=all(checks.values())
print('EXACT PLANE CHECKS',len(checks),all(checks.values()),flush=True)
plane_vector_fields=[[pull(src['Fraw'][i][j]).subs(exact_plane) for j in [2,5]] for i in range(3)]
summary['plane_fields_pq']=[[str(z) for z in row] for row in plane_vector_fields]
print('PLANE FIELDS',summary['plane_fields_pq'],flush=True)
save({'ring':R,'leading_pi':ps,'t_graph':sub,'first_residuals':res,'equations':eq},str(OUT/'homothety_axis_infinity_leading.sobj'))
(OUT/'homothety_axis_infinity_leading.json').write_text(json.dumps(summary,indent=2))
print('COMPLETE',round(time.time()-START,3),flush=True)
