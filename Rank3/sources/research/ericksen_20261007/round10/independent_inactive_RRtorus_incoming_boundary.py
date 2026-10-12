"""Literal full-source first boundary on the inactive full-RR incoming support.

Only a source inventory and regular finite-range chart. No no-lift assertion.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(s,b):assert b,s;checks.append(s)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset'];hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'));sp=load(str(p/'independent_supercritical_first_boundary_probe.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
K=qf['ring'].base_ring();a,r,H,Z=K.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));shift=-2*H*(Z+1)/(Z*(H+1));P=PolynomialRing(K,names=('in_l','in_rho','in_e0','in_e1','in_e2','in_delta'));l,rho,e0,e1,e2,delta=P.gens();PF=P.fraction_field();es=[e0,e1,e2];vals=[P.one(),0,0,0,0,l,alpha*l,rho,P.one()]
def pull_dict(dd,cv=lambda x:K(x)):
 return sum((P(cv(co))*prod(P(vals[j])**int(ex[j]) for j in range(9)) for ex,co in dd.items() if not any(vals[j]==0 and ex[j] for j in range(9))),P.zero())
def pull(f):return pull_dict(f.dict())
def layer(f,bi):return {tuple(map(int,e)):co for e,co in f.dict().items() if (sum(e[:6]),sum(e[6:]))==bi}
C=matrix(P,3,3,[pull(f) for f in qf['C'].list()]);D=matrix(P,6,3,[pull(f) for f in qf['D'].list()]);dd=vector(P,[pull(f) for f in off]);Rrows=list(D*(C*vector(P,es)+dd)+delta*vector(P,[pull(f) for f in sp['F12']]));Qrow=sum((pull(f)*es[j] for j,f in enumerate(qf['Q0_L_cofactors'])),P.zero())+pull_dict(layer(qf['full_Q0'],(2,1)));rows=Rrows+[Qrow];records=[('R',j) for j in range(6)]+[('Q0',0)]
print('R Q0 formed',time.monotonic()-st,flush=True)
old=load(str(p/'independent_delta_H_global_L_module.sobj'))['original_field'];km=old.ring().hom([a,a*(1+r)/(1+a*r),H/r,(1+r)*Z/r],K);cache={}
def cv(z):
 z=old(z)
 if z not in cache:cache[z]=km(z.numerator())/km(z.denominator())
 return cache[z]
for j,it in enumerate(hi['items']):
 ans=sum((pull_dict(dd,cv)*es[k] for k,dd in enumerate(it['L_cofactors_original'])),P.zero())+pull_dict(it['H41_original'],cv)+delta*pull(sp['H13'][j]);rows.append(ans);records.append(('H',j));print('H formed',j,time.monotonic()-st,flush=True)
even=[]
for it in ev['items']:
 even.append(K(it['square_unit'])*es[int(it['square_component'])]**2+sum((pull(f)*es[j] for j,f in enumerate(it['A_p_quadratic'])),P.zero())+pull(it['pure_p4']))
sub={e0:shift*l,e1:0,e2:0,delta:0};res=[f.subs(sub) for f in rows+even]
for j,f in enumerate(res):ck('whole_incoming_corner_zero_'+str(j),f==0)
M=matrix(P,[[f.derivative(e0),f.derivative(e1),f.derivative(e2),f.derivative(delta),f.subs({e0:0,e1:0,e2:0,delta:0})] for f in rows]);ck('affine_matrix_reconstruct',M*vector(P,[e0,e1,e2,delta,1])==vector(P,rows))
print('M shape',M.dimensions(),flush=True)
print('EVEN term counts',[len(f.dict()) for f in even],flush=True)
print('RROWS term counts',[len(f.dict()) for f in rows[:7]],flush=True)
# Chart identity, preserving both ratio parameters and all three range variables.
U=PolynomialRing(K,names=('at_l','at_rho','at_s','at_w','at_e0','at_e1','at_e2','at_k','at_sigma'));ll,rr,ss,ww,E0,E1,E2,k,sig=U.gens();UF=U.fraction_field();vv=[sig*(alpha*ll+k*E0),sig*rr,sig];pp=[1,((alpha/beta)*ss*rr+k*E1)/(alpha*ll+k*E0),(beta*ww+k*E2)/rr,ss,ww,ll];tau=k*sig
ells=[pp[0]*vv[0]-alpha*pp[5]*vv[2],pp[1]*vv[0]-(alpha/beta)*pp[3]*vv[1],pp[2]*vv[1]-beta*pp[4]*vv[2]]
for j,f in enumerate(ells):ck('exact_finite_range_chart_'+str(j),UF(f)==tau*[E0,E1,E2][j])
ck('incoming_ordinary_relation',UF(vv[0]/tau).subs({E0:shift*ll})==ll*(alpha/k+shift))
save({'coefficient_field':K,'ring':P,'p_RR_images':vals,'shift':shift,'C':C,'D':D,'rows':rows,'records':records,'matrix':M,'even':even,'incoming_range_graph':sub,'chart_ring':U,'chart_p_images':pp,'chart_RR_images':vv,'chart_tau':tau,'checks':checks,'scope':__doc__},str(p/'independent_inactive_RRtorus_incoming_boundary.sobj'),compress=False)
(p/'independent_inactive_RRtorus_incoming_boundary.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
