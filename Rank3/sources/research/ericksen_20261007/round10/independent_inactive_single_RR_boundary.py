"""Full original intermediate-boundary equations on all three single-RR axes.

Retain all four allowed p amplitudes and all finite-range variables.
This script constructs necessary equations, without a no-lift assertion.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):assert v,n;checks.append(n)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset']
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
sp=load(str(p/'independent_supercritical_first_boundary_probe.sobj'))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
K=qf['ring'].base_ring();a,r,H,Z=K.gens()
P=PolynomialRing(K,names=('ax_x','ax_y','ax_z','ax_w','ax_e0','ax_e1','ax_e2','ax_delta'))
x,y,z,w,e0,e1,e2,de=P.gens();es=[e0,e1,e2]
old=load(str(p/'independent_delta_H_global_L_module.sobj'))['original_field']
km=old.ring().hom([a,a*(1+r)/(1+a*r),H/r,(1+r)*Z/r],K);cache={}
def cv(c):
 c=old(c)
 if c not in cache:cache[c]=km(c.numerator())/km(c.denominator())
 return cache[c]
def layer(f,bi):return {tuple(map(int,e)):co for e,co in f.dict().items() if (sum(e[:6]),sum(e[6:]))==bi}
cases=[(0,0,x,y,z,w,1,0,0),(x,y,0,0,z,w,0,1,0),(x,y,z,w,0,0,0,0,1)]
items=[]
for ai,vals in enumerate(cases):
 def pd(dd,coerce=lambda c:K(c)):
  out=P.zero()
  for ex,c in dd.items():
   if any(vals[j]==0 and ex[j] for j in range(9)):continue
   out+=P(coerce(c))*prod(P(vals[j])**int(ex[j]) for j in range(9))
  return out
 def pull(f):return pd(f.dict())
 C=matrix(P,3,3,[pull(f) for f in qf['C'].list()]);D=matrix(P,6,3,[pull(f) for f in qf['D'].list()]);dd=vector(P,[pull(f) for f in off])
 Rrows=list(D*(C*vector(P,es)+dd)+de*vector(P,[pull(f) for f in sp['F12']]))
 Qrow=sum((pull(f)*es[j] for j,f in enumerate(qf['Q0_L_cofactors'])),P.zero())+pd(layer(qf['full_Q0'],(2,1)))
 rows=Rrows+[Qrow]
 for j,it in enumerate(hi['items']):
  rows.append(sum((pd(f,cv)*es[k] for k,f in enumerate(it['L_cofactors_original'])),P.zero())+pd(it['H41_original'],cv)+de*pull(sp['H13'][j]))
 even=[]
 for it in ev['items']:
  even.append(K(it['square_unit'])*es[int(it['square_component'])]**2+sum((pull(f)*es[j] for j,f in enumerate(it['A_p_quadratic'])),P.zero())+pull(it['pure_p4']))
 M=matrix(P,[[f.derivative(e0),f.derivative(e1),f.derivative(e2),f.derivative(de),f.subs({e0:0,e1:0,e2:0,de:0})] for f in rows])
 ck('affine_reconstruct_'+str(ai),M*vector(P,[e0,e1,e2,de,1])==vector(P,rows))
 print('AXIS',ai,'D',D,'R Q0',*[str(f.factor()) if f else '0' for f in rows[:7]],'EVEN',*[str(f.factor()) if f else '0' for f in even],sep='\n',flush=True)
 print('H sizes',[len(f.dict()) for f in rows[7:]],'seconds',time.monotonic()-st,flush=True)
 items.append({'axis':ai,'images':vals,'D':D,'C':C,'rows':rows,'matrix':M,'even':even})
 save({'ring':P,'items':items,'checks':checks,'scope':__doc__},str(p/'independent_inactive_single_RR_boundary.sobj'),compress=False)
print('DONE',len(checks),'seconds',time.monotonic()-st,flush=True)
