"""Exact ingredients for weighted higher-jet bounds near a principal axis."""
from sage.all import *
from pathlib import Path
import pickle,json,time
OUT=Path(__file__).absolute().parent;START=time.time()
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
K=FractionField(PolynomialRing(QQ,names=('a','b','s','t')));a,b,s,t=K.gens()
R=PolynomialRing(K,names=tuple('p'+str(i) for i in range(9))+tuple('v'+str(i) for i in range(6))+('u0','u1'))
ps=R.gens()[:9];vs=R.gens()[9:15];us=R.gens()[15:];n=(s,t,1)
def import_poly(z): return R({ex:K(cf) for ex,cf in z.items()})
F=[[import_poly(z) for z in row] for row in src['Fraw']]
H=[sum(n[i]*F[i][j] for i in range(3)) for j in range(17)]
mt=matrix(K,[[z.monomial_coefficient(v) for v in vs] for z in H[:9]])
mu=matrix(K,[[z.monomial_coefficient(u) for u in us] for z in H[9:15]])
mt0=mt.apply_map(lambda z: z.subs({s:0,t:0}));mu0=mu.apply_map(lambda z:z.subs({s:0,t:0}))
assert mt0.rank()==5 and mu0.rank()==2
tr=list(mt0.transpose().pivots());tc=list(mt0.pivots());ur=list(mu0.transpose().pivots())
assert tc==[1,2,3,4,5]
assert mu.matrix_from_rows(ur)==matrix.identity(K,2)
tdet=mt0.matrix_from_rows_and_columns(tr,tc).det();udet=mu0.matrix_from_rows(ur).det()
zpi={p:0 for p in ps};zu={u:0 for u in us};zv={v:0 for v in vs[1:]}
q4=[import_poly(z) for z in src['Q4']]
pure_t=[q.subs({**zpi,**zu,**zv}) for q in q4]
indices=[i for i,q in enumerate(pure_t) if q]
qindex=indices[0]; coefficient=pure_t[qindex].monomial_coefficient(vs[0]**2)
assert coefficient==8 and pure_t[qindex]==8*vs[0]**2
weights=[1]*9+[2]*6+[3]*2
for i,row in enumerate(F):
    for j,q in enumerate(row):
        assert all(sum(ex[k]*weights[k] for k in range(17))==weights[j]+1 for ex in q.dict())
assert all(sum(ex[k]*weights[k] for k in range(17))==4 for ex in q4[qindex].dict())
# Weight four prohibits u*u or u*t, and makes every u term carry a pi.
assert all(sum(ex[15:])<=1 and (not sum(ex[15:]) or sum(ex[:9])==1) for ex in q4[qindex].dict())
rec={'status':'EXACT_INGREDIENTS_FOR_WEIGHTED_JET_BOUND', 'axis_n':[0,0,1],
 't_rows':tr,'t_columns':tc,'t_pivot':str(tdet.factor()),
 'u_rows':ur,'u_pivot':str(udet.factor()),'u_block_identity_for_every_direction':True,'u_matrix_at_axis':[[str(z) for z in row] for row in mu0.rows()],
 'Q4_index':qindex,'Q4_on_pi0_u0_t_except0zero':str(pure_t[qindex]),
 'all_checks_pass':True,'seconds':time.time()-START}
print(json.dumps(rec,indent=2),flush=True)
save({'ring':R,'t_matrix':mt,'u_matrix':mu,'t_rows':tr,'t_columns':tc,'u_rows':ur,'quartic':q4[qindex],'Q4_index':qindex},str(OUT/'homothety_axis_jet_bound_reader.sobj'))
(OUT/'homothety_axis_jet_bound_reader.json').write_text(json.dumps(rec,indent=2))
