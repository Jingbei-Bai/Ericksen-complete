"""Original-source normal linearization along the entire axis infinity plane.

No graph with 1/s or 1/t is used. This is a necessary tangent calculation,
not exclusion of incoming arcs. Parameters p,q are not specialized.
"""
from sage.all import *
from pathlib import Path
import json,pickle,time
OUT=Path(__file__).absolute().parent;START=time.time()
src=pickle.loads((OUT.parent/'round3/generic_degree5_sparse_prepared.sobj').read_bytes())
K=FractionField(PolynomialRing(QQ,names=('a','b','p','q')));a,b,p,q=K.gens()
names=tuple('r'+str(i) for i in range(9))+tuple('v'+str(i) for i in range(6))+('u0','u1')
R=PolynomialRing(K,names=names);x=R.gens();weights=[1]*9+[2]*6+[3]*2
base=vector(K,[0,0,p,0,0,q,0,0,0,p*q,0,0,0,0,0,0,0])
F=[[R({ex:K(cf) for ex,cf in comp.items()}) for comp in row] for row in src['Fraw']]
def original(z):return R({ex:K(cf) for ex,cf in z.items()})
def at(f):return K(f(*base))
S=PolynomialRing(K,names=tuple('d'+str(i) for i in range(17))+('epsilon','n1','n2'))
dx=S.gens()[:17];eps,n1,n2=S.gens()[17:]
normal=[];labels=[];base_checks={}
for j,z in enumerate(F[2]):
    base_checks[f'H_{j}']=at(z)==0
    ll=sum(at(z.derivative(x[k]))*dx[k] for k in range(17))
    ll+=weights[j]*base[j]*eps+at(F[0][j])*n1+at(F[1][j])*n2
    normal.append(ll);labels.append(('H',j))
for name,seq in [('C3',src['C3']),('Q4',src['Q4']),('G5',src['basis'][:27])]:
    for j,comp in enumerate(seq):
        z=original(comp);base_checks[f'{name}_{j}']=at(z)==0
        ll=sum(at(z.derivative(x[k]))*dx[k] for k in range(17))
        normal.append(ll);labels.append((name,j))
    print('normal',name,len(normal),'seconds',time.time()-START,flush=True)
assert all(base_checks.values()),[k for k,v in base_checks.items() if not v]
M=matrix(K,[[z.monomial_coefficient(y) for y in S.gens()] for z in normal])
N=M[:,:17];rr=list(N.transpose().pivots());cc=list(N.pivots());minor=N.matrix_from_rows_and_columns(rr,cc)
ker=N.right_kernel().basis_matrix()
pi=matrix(R,3,3,x[:9]); lam=diagonal_matrix(K,[1,1/a,1/b]); I=matrix.identity(R,3)
nbase=vector(R,[0,0,1]); pressure=matrix(R,[[pi[1,2]-pi[2,1],pi[0,2],-pi[0,1]],
 [-pi[1,2],pi[2,0]-pi[0,2],pi[1,0]],[pi[2,1],-pi[2,0],pi[0,1]-pi[1,0]]])
obs=list((pi.transpose()-pi.trace()*I/2)*nbase)+list((pi.transpose()*lam-(lam*pi).trace()*I/2)*nbase)+list((pressure*nbase)[:2])
O=matrix(K,[[at(z.derivative(x[k])) for k in range(17)] for z in obs])
stack=O.stack(N); stackrows=list(stack.transpose().pivots());stackcols=list(stack.pivots())
stackminor=stack.matrix_from_rows_and_columns(stackrows,stackcols).det()
strata={}
for label,substitution in [('p_nonzero_q_zero',{q:0}),('p_zero_q_nonzero',{p:0})]:
 NN=N.apply_map(lambda z:z.subs(substitution));SS=O.apply_map(lambda z:z.subs(substitution)).stack(NN)
 r1=list(NN.transpose().pivots());c1=list(NN.pivots());r2=list(SS.transpose().pivots());c2=list(SS.pivots())
 strata[label]={'source_rank':len(c1),'observations_and_source_rank':len(c2),
  'source_minor':str(NN.matrix_from_rows_and_columns(r1,c1).det().factor()),
  'combined_minor':str(SS.matrix_from_rows_and_columns(r2,c2).det().factor())}
print('normal rank',len(cc),'cols',cc,'rows',rr,'minor',minor.det().factor(),flush=True)
print('kernel',ker,flush=True)
summary={'status':'NECESSARY_LINEARIZATION_ONLY','base_checks':base_checks,'normal_rank':len(cc),
 'rows':rr,'columns':cc,'pivot_determinant':str(minor.det().factor()),
 'kernel_rows':[[str(z.factor()) if z else '0' for z in row] for row in ker.rows()],
 'normal_counts':{name:sum(z!=0 for z,lab in zip(normal,labels) if lab[0]==name) for name in ('H','C3','Q4','G5')},
 'observations_plus_normal_rank':len(stackcols),'combined_minor':str(stackminor.factor()),'boundary_strata':strata,
 'seconds':float(time.time()-START)}
save({'ring':S,'parameter_field':K,'base_state':base,'weights':weights,'normal_sources':normal,'labels':labels,
      'matrix':M,'state_matrix':N,'kernel':ker,'pivot_rows':rr,'pivot_columns':cc,'source_fields':F,
      'source_ring':R,'base_checks':base_checks,'observation_matrix':O,
      'observation_plus_normal_matrix':stack,'combined_kernel':stack.right_kernel().basis_matrix(),
      'combined_pivot_rows':stackrows,'combined_pivot_columns':stackcols,
      'combined_pivot_determinant':stackminor},str(OUT/'homothety_axis_plane_normal_sources.sobj'))
(OUT/'homothety_axis_plane_normal_sources.json').write_text(json.dumps(summary,indent=2))
print('COMPLETE',summary['seconds'],flush=True)
