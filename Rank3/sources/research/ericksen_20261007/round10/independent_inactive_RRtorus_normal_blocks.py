"""Small genuine block minors at the inactive incoming RR-torus infinity."""
from sage.all import *
from pathlib import Path
from itertools import combinations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_inactive_RRtorus_incoming_boundary.sobj'));P=d['ring'];l,rho,e0,e1,e2,de=P.gens();K=P.base_ring();a,r,H,Z=K.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));shift=d['shift'];M=d['matrix'];Q=PolynomialRing(K,names=('bl','brho'));L,R=Q.gens();QF=Q.fraction_field()
def conv(f):
 f=P(f);return Q({tuple(map(int,ex[:2])):co for ex,co in f.dict().items()})
N=matrix(Q,25,3,[conv(M[i,j]) for i in range(25) for j in [1,2,3]])
seed=N.matrix_from_rows_and_columns([1,3],[0,1]);print('seed nonzero',bool(seed.det()),flush=True)
rs=[]
for j in range(7,25):
 if not N.row(j):continue
 v=N.matrix_from_rows([1,3,j]).det();print('DET',j,'terms',len(v.dict()),flush=True);rs.append({'row':j,'minor':v})
q=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));O=q['ring'];xs=O.gens();ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));vals=[Q.one(),0,0,0,0,L,0,0,0]
def pull(f):return sum((Q(co)*prod(Q(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items() if not any(vals[j]==0 and ex[j] for j in range(9))),Q.zero())
evrows=[]
for it in ev['items']:
 # Derivatives along paired-support transverse p coordinates, with e fixed at its true graph.
 dp=[]
 for ix in [1,3,2,4]:
  v=pull(it['pure_p4'].derivative(xs[ix]))
  v+=pull(it['A_p_quadratic'][0].derivative(xs[ix]))*shift*L
  dp.append(v)
 evrows.append([QF(R/(beta*L))*dp[0]+dp[1],QF(beta/R)*dp[2]+dp[3]])
E=matrix(QF,evrows);emin=[]
for ij in combinations(range(3),2):
 vv=E.matrix_from_rows(ij).det();print('EVEN_DET',ij,'nonzero',bool(vv),flush=True);emin.append({'rows':ij,'minor':vv})
shortP=R**2/beta**2-5*r*L**2-4*(1+r);shortQ=R**2/beta**2+4*r*L**2+5*(1+r);m=1+a*r
vp=next(z['minor'] for z in rs if z['row']==8);vq=next(z['minor'] for z in rs if z['row']==20)
up=QF(vp)/(L**4*shortP);uq=QF(vq)/(L**3*shortQ);ue=emin[0]['minor']*R/(L*((m*L**2-a*r)*(m*L**2+a*r)))
for z0 in [up,uq,ue]:assert z0.numerator().total_degree()==0 and z0.denominator().total_degree()==0
assert shortQ-shortP==9*(r*L**2+1+r)
short={'P':shortP,'Q':shortQ,'P_source_row':8,'Q_source_row':20,'P_unit':up,'Q_unit':uq,'even_unit':ue,'even_factor':(m*L**2-a*r)*(m*L**2+a*r),'Q0_e0':conv(M[6,0])}
print('UNITS',[str(K(z0).factor()) for z0 in [up,uq,ue]],flush=True)
save({'short_certificates':short,'source_file':'independent_inactive_RRtorus_incoming_boundary.sobj','ring':Q,'range_delta_block':N,'seed_rows':[1,3],'seed_minor':seed.det(),'H_minors':rs,'even_transverse_matrix':E,'even_minors':emin,'scope':__doc__},str(p/'independent_inactive_RRtorus_normal_blocks.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
