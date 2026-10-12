"""Entire nonzero-scale common-homothety first-jet reduction.

Original 17-state sources, all material parameters. sigma=1 is a common
generator normalization; u3!=0 is justified by the principal-pressure identity.
No zero of this finite necessary system is asserted to be an actual solution.
"""
from sage.all import *
from pathlib import Path
import pickle,time,json,resource
OUT=Path(__file__).absolute().parent;OUT.mkdir(exist_ok=True);BASE=OUT.parent;START=time.time()
source=pickle.loads((BASE/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
BR=PolynomialRing(QQ,names=('a','b','s','t'));a,b,s,t=BR.gens();K=BR.fraction_field()
names=['p11','p12','p13','p21','p22','p23','p31','p32','p33']+[f't{i}' for i in range(6)]+['u0','u1']
R=PolynomialRing(K,names=names);vars=R.gens();P=vars[:9];ts=vars[9:15]
F=[[R({ex:K(cf) for ex,cf in comp.items()}) for comp in row] for row in source['Fraw']]
pi=matrix(R,3,3,P);n=vector(K,[s,t,1])
eq=vector(R,[sum(n[i]*F[i][q] for i in range(3))+P[q] for q in range(9)])
M=matrix(K,[[z.monomial_coefficient(tt) for tt in ts] for z in eq]);zsub={tt:0 for tt in ts}
rhs=vector(R,[-z.subs(zsub) for z in eq])
assert all(z.derivative(tt)==M[i,j] for i,z in enumerate(eq) for j,tt in enumerate(ts))
rows=[0,1,2,3,5,6];cols=list(M.pivots());block=M.matrix_from_rows_and_columns(rows,cols)
print('rank',len(cols),'rows',rows,'cols',cols,'det',block.det().factor(),flush=True)
if len(cols)!=6:raise RuntimeError('first homothety derivatives do not determine all six t')
sol=block.inverse()*vector(R,[rhs[i] for i in rows]);res=M*sol-rhs
for ix,z in enumerate(res):
 if z:print('res',ix,'terms',len(z.dict()),'degree',z.total_degree(),flush=True)
def eps(i,j,k):
 if len({i,j,k})<3:return 0
 return 1 if (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)] else -1
c={(l,k,i):-sum(eps(k,i,m)*pi[l,m] for m in range(3)) for l in range(3) for k in range(3) for i in range(3)}
gamma=vector(R,[int(i==2)+sum(c[2,k,i]*n[k] for k in range(3)) for i in range(3)])
dn=[vector(R,[int(i==h)+sum(c[h,k,i]*n[k] for k in range(3))-n[h]*gamma[i] for h in range(2)]) for i in range(3)]
fp=[[z.subs(dict(zip(ts,sol)))+gamma[i]*P[h] for h,z in enumerate(F[i][:9])] for i in range(3)]
assert sum(n[i]*gamma[i] for i in range(3))==1
pressure=matrix(R,[[pi[1,2]-pi[2,1],pi[0,2],-pi[0,1]],[-pi[1,2],pi[2,0]-pi[0,2],pi[1,0]],[pi[2,1],-pi[2,0],pi[0,1]-pi[1,0]]])
source_omega=(pi.transpose()-pi.trace()/2*matrix.identity(R,3))*vector(R,n)
lam=diagonal_matrix(K,[1,1/a,1/b]);tlam=2*pi.transpose()*lam-(lam*pi).trace()*matrix.identity(R,3)
target_omega_norm=R((vector(R,n)*tlam.transpose()*lam*tlam*vector(R,n))[0]/(4*lam.det())) if False else sum(lam[i,i]*(tlam*vector(R,n))[i]**2 for i in range(3))/(4*lam.det())
save({'ring':R,'parameters':(a,b,s,t),'pi':pi,'direction':n,'fields':F,'matrix':M,'rows':rows,'determinant':block.det(),'t_solution':sol,'residuals':res,'gamma':gamma,'direction_fields':dn,'scaled_pi_fields':fp,'longitudinal_pressures':pressure*vector(R,n),'source_rotation_vector':source_omega,'source_rotation_norm':source_omega.dot_product(source_omega),'target_rotation_norm':target_omega_norm},str(OUT/'homothety_first_jet_reduction.sobj'))
summary={'status':'EXACT_NECESSARY_REDUCTION','rows':rows,'determinant':str(block.det().factor()),'residuals':[(i,str(z.factor())) for i,z in enumerate(res) if z],'seconds':time.time()-START}
(OUT/'homothety_first_jet_reduction.json').write_text(json.dumps(summary,indent=2))
print('COMPLETE',time.time()-START,flush=True)
