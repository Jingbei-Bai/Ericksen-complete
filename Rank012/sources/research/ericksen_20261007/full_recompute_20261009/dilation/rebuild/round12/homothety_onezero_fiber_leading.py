"""All-material first-jet fiber leaders on a single zero direction wall."""
from sage.all import *
from pathlib import Path
import pickle,json,time
OUT=Path(__file__).absolute().parent; START=time.time()
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
K=FractionField(PolynomialRing(QQ,names=('a','b','v')));a,b,v=K.gens()
R=PolynomialRing(K,names=tuple('p'+str(i) for i in range(9))+tuple('t'+str(i) for i in range(6))+('u0','u1'))
ps=R.gens()[:9];ts=R.gens()[9:15];pi=matrix(R,3,3,ps)
F=[[R({ex:K(cf) for ex,cf in comp.items()}) for comp in row] for row in src['Fraw']]
Lam=diagonal_matrix(K,[1,1/a,1/b]);Id=matrix.identity(R,3)
P=matrix(R,[[pi[1,2]-pi[2,1],pi[0,2],-pi[0,1]],[-pi[1,2],pi[2,0]-pi[0,2],pi[1,0]],[pi[2,1],-pi[2,0],pi[0,1]-pi[1,0]]])
summary=[];records=[]
for direction in [[0,v,1],[v,0,1]]:
    n=vector(K,direction)
    obs=list((pi.transpose()-pi.trace()/2*Id)*n)+list((pi.transpose()*Lam-(Lam*pi).trace()/2*Id)*n)+list((P*n)[:2])
    M=matrix(K,[[z.monomial_coefficient(p) for p in ps] for z in obs]); assert M.rank()==8
    kernel=M.right_kernel().basis()[0]
    assert kernel[0]
    kernel=kernel/kernel[0]
    eq=vector(R,[sum(n[i]*F[i][j] for i in range(3)) for j in range(9)])
    T=matrix(K,[[z.monomial_coefficient(tt) for tt in ts] for z in eq]);assert T.rank()==6
    rows=list(T.transpose().pivots());block=T.matrix_from_rows(rows)
    pivot=block.det()
    # First make the fiber kernel have amplitude p0; retain p0 as a free symbol.
    zsub={ps[j]:R(kernel[j])*ps[0] for j in range(9)}
    rhoeq=vector(R,[z.subs(zsub) for z in eq])
    rhs=vector(R,[-rhoeq[i].subs({tt:0 for tt in ts}) for i in rows])
    sol=block.inverse()*rhs
    remaining=[z.subs(dict(zip(ts,sol))) for z in rhoeq]
    leaders=[K(z.monomial_coefficient(ps[0]**2)) for z in remaining if z]
    assert all(z==z.monomial_coefficient(ps[0]**2)*ps[0]**2 for z in remaining)
    rec={'direction':[str(z) for z in n], 'kernel':[str(z.factor()) if z else '0' for z in kernel],
       't_pivot':str(pivot.factor()), 'leaders':[str(z.factor()) for z in leaders]}
    print(json.dumps(rec,indent=2),flush=True);summary.append(rec)
    records.append({'direction':n,'observation_matrix':M,'kernel':kernel,'t_matrix':T,'t_rows':rows,'t_pivot':pivot,'leaders':leaders})
save({'ring':R,'records':records},str(OUT/'homothety_onezero_fiber_leading.sobj'))
(OUT/'homothety_onezero_fiber_leading.json').write_text(json.dumps({'status':'NECESSARY_FIBER_LEADERS_NOT_CLASSIFICATION','records':summary,'seconds':time.time()-START},indent=2))
print('COMPLETE',round(time.time()-START,3),flush=True)
