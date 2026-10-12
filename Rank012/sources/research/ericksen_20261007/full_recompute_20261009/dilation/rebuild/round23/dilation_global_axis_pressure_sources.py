"""Bare principal-axis necessary source; no selected-spectrum positive inference."""
from sage.all import *
from pathlib import Path
import pickle,json,time,signal
signal.alarm(60)
O=Path(__file__).absolute().parent;start=time.time()
src=pickle.loads((O.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
K=FractionField(PolynomialRing(QQ,names=('a','b')));a,b=K.gens()
names=('h1','p12','p13','p21','h2','p23','p31','p32','h3')+tuple('t'+str(i) for i in range(6))+('u0','u1')
R=PolynomialRing(K,names=names);y=R.gens();ps=y[:9];ts=y[9:15];us=y[15:]
def pull(comp):
 return sum(R(K(cf))*prod(y[i]**ex[i] for i in range(17)) for ex,cf in comp.items())
F=[[pull(c) for c in row] for row in src['Fraw']]
hp=vector(R,[F[2][i]+ps[i] for i in range(9)])
M=matrix(K,[[z.monomial_coefficient(tt) for tt in ts] for z in hp])
ri=list(M.transpose().pivots());ci=list(M.pivots())
rhs=vector(R,[-hp[i].subs({ts[j]:0 for j in ci}) for i in ri])
sol=M.matrix_from_rows_and_columns(ri,ci).inverse()*rhs
sub={ts[j]:sol[i] for i,j in enumerate(ci)}
res=[z.subs(sub) for z in hp if z.subs(sub)]
ht=[(F[2][i]+2*y[i]).subs(sub) for i in range(9,15)]
c3=[pull(c).subs(sub) for c in src['C3']]
rep={'status':'NECESSARY_AXIS_SOURCES_ONLY','rank':M.rank(),'pivotrows':ri,'pivotcols':ci,
 'pivotdet':str(M.matrix_from_rows_and_columns(ri,ci).det().factor()),
 'Hpi_residuals':[str(z.factor()) for z in res],
 'Hpi_terms':[len(z.dict()) for z in res],
 'Ht_terms':[len(z.dict()) for z in ht], 'C3_terms':[len(z.dict()) for z in c3],
 'seconds':time.time()-start}
save({'ring':R,'fields':F,'Hpi':hp,'t_graph':sub,'Hpi_residuals':res,'Ht':ht,'C3':c3},str(O/'dilation_global_axis_pressure_sources.sobj'))
(O/'dilation_global_axis_pressure_sources.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
print(json.dumps(rep,indent=2),flush=True)
