"""Exact complete regular sources on the p13=1 single RR incidence chart.

No truncated derivatives. Replace only certified leading L modules;
every lower original coefficient, h/z derivative, and positive Newton
power is retained. Specialization of coefficients is an optional probe.
"""
from sage.all import *
from pathlib import Path
import time,sys
p=Path(__file__).resolve().parent;st=time.monotonic();qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));R=qf['ring'];OK=R.base_ring();oa,orr,oH,oZ=OK.gens();fixed=len(sys.argv)>1
K=PolynomialRing(QQ,names=('fc_k','fc_r') if fixed else ('fc_k','fc_r','fc_H','fc_Z')).fraction_field();k,r=K.gens()[:2];H,Z=(K(2),K(-1)) if fixed else K.gens()[2:];m=1+k;q=1+r;ph=OK.ring().hom([k/r,r,H,Z],K);cache={}
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
P=PolynomialRing(K,names=('fc_x','fc_z','fc_w','fc_e0','fc_e1','fc_e2','fc_delta','fc_t','fc_kappa'));x,z,w,e0,e1,e2,de,t,kap=P.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));imgs=(x,P.one(),z,w,(t*z-kap*e2)/beta,x*w*t/beta+kap*(x*e1-e0)/alpha,alpha/beta*w*t+kap*e1,t,P.one());es=(e0,e1,e2)
def cp(f):return R({tuple(map(int,ee))+(0,)*(9-len(ee)):OK(c) for ee,c in f.dict().items()})
def mp(f):return sum((cv(c)*prod(v**int(n) for v,n in zip(imgs,ee)) for ee,c in cp(f).dict().items()),P.zero())
def norm(f,offset,discard):
 ans=P.zero()
 for ee,c in cp(f).dict().items():
  pd=sum(ee[:6]);rd=sum(ee[6:]);
  if (pd,rd) in discard:continue
  ke,dd=offset[0]-2*pd-3*rd,offset[1]-pd-rd
  assert ke>=0 and dd>=0,(pd,rd,ke,dd)
  ans+=cv(c)*kap**ke*de**dd*prod(v**int(n) for v,n in zip(imgs,ee))
 return ans
off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset'];off=vector(R,[cp(f) for f in off]);DC=matrix(R,qf['D'])*matrix(R,qf['C']);sources=[];labels=[]
for j,full in enumerate(qf['full_R_sources']):sources.append(norm(full,(8,4),{(3,1)})+sum((mp(DC[j,i])*es[i] for i in range(3)),P.zero()));labels.append(('R',j))
sources.append(norm(qf['full_Q0'],(7,3),{(1,2)})+sum((mp(qf['Q0_L_cofactors'][i])*es[i] for i in range(3)),P.zero()));labels.append(('Q0',0));print('R/Q built',time.monotonic()-st,flush=True)
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));old=load(str(p/'independent_delta_H_global_L_module.sobj'))['original_field'];oph=old.ring().hom([k/r,(k/r)*q/m,H/r,q*Z/r],K);oldcache={};pool=full['polynomial_pool']
def oldcv(c):
 c=old(c)
 if c not in oldcache:oldcache[c]=oph(c.numerator())/oph(c.denominator())
 return oldcache[c]
def om(f):return sum((oldcv(c)*prod(v**int(n) for v,n in zip(imgs,ee)) for ee,c in f.items()),P.zero())
for j,it in enumerate(hs['items']):
 coeffs={}
 for sign,lr,rr in full['records'][j]['full_signed_product_DAG']:
  for le,lc in pool[lr].items():
   for re,rc in pool[rr].items():
    ee=tuple(int(v)+int(u) for v,u in zip(le,re))
    if (sum(ee[:6]),sum(ee[6:]))==(3,2):continue
    coeffs[ee]=coeffs.get(ee,K.zero())+sign*oldcv(lc)*oldcv(rc)
 coeffs={ee:cc for ee,cc in coeffs.items() if cc};ans=sum((om(it['L_cofactors_original'][i])*es[i] for i in range(3)),P.zero())
 for ee,c in coeffs.items():
  pd=sum(ee[:6]);rd=sum(ee[6:]);
  if (pd,rd)==(3,2):continue
  ke,dd=11-2*pd-3*rd,5-pd-rd;assert ke>=0 and dd>=0,(j,pd,rd,ke,dd)
  ans+=c*kap**ke*de**dd*prod(v**int(n) for v,n in zip(imgs,ee))
 sources.append(ans);labels.append(('H',j));print('H',j,'terms',len(ans.dict()),'sec',time.monotonic()-st,flush=True)
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for it in ev['items']:
 src=cp(it['full_source']);ell=matrix(R,qf['L'])*vector(R,R.gens()[6:]);rem=src-OK(it['square_unit'])*ell[int(it['square_component'])]**2-sum((cp(it['A_p_quadratic'][i])*ell[i] for i in range(3)),R.zero())
 ans=norm(rem,(8,4),set())+cv(it['square_unit'])*es[int(it['square_component'])]**2+sum((mp(it['A_p_quadratic'][i])*es[i] for i in range(3)),P.zero());sources.append(ans);labels.append(('even',it['Q_index']))
bd=load(str(p/'independent_critical_single_RR3_normal.sobj'));BP=bd['ring'];zero=P.hom([x,z,w,e0,e1,e2,de,P.zero(),P.zero()],P)
def bm(f):return sum((cv(c)*prod(v**int(n) for v,n in zip(P.gens(),ee)) for ee,c in f.dict().items()),P.zero())
for j,f in enumerate(sources):assert zero(f)==bm(bd['boundary'][j]),j
stem='independent_subcritical_single_RR_full_regular_chart'+('_H2' if fixed else '')
save({'ring':P,'sources':sources,'labels':labels,'chart_p':imgs[:6],'chart_Rbar':imgs[6:],'source_offset_exponents':{'R':(8,4),'Q0':(7,3),'H':(11,5),'even':(8,4)},'scope':__doc__},str(p/(stem+'.sobj')),compress=False);print('PASS whole boundary',len(sources),'sec',time.monotonic()-st,flush=True)
