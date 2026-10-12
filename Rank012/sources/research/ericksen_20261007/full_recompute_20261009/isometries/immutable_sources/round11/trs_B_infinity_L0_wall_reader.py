"""Whole local-ring all-pitch certificate at the exact determinant wall.

Factor guards give regular sources; source-linked jets check the finite
quotients (x,tau)/(x,tau)^2 and J/mJ, J=(x,tau^2).
"""
from sage.all import *
from pathlib import Path
import json,resource
OUT=Path(__file__).resolve().parent
d=load(str(OUT/'trs_algebra_regular_frobenius.sobj'));F=d['source']['field'];PR=F.ring()
v,p,j,A,T,B,L=map(F,PR.gens())
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
checks={};H=max((ff for ff,ee in d['frobenius'][1][0][(0,0)].denominator().factor()),key=lambda ff:len(ff.dict()))
checks['wall_polynomial_linear_L']=H.degree(PR.gen(6))==1
allowed=[PR.gen(k) for k in [0,2,3,4,5]]+[PR.gen(0)+PR.gen(6)**2,PR.gen(0)-PR.gen(1)**2]
for row,C in enumerate(d['frobenius']):
 for odd,pol in enumerate(C):
  for press,f in pol.items():
   if not f:continue
   rem=F(f.denominator())/F(H)**3;ok=True
   for poly in [rem.numerator(),rem.denominator()]:
    for fac,mult in poly.factor():
     # Sage associates may carry a sign only.
     if not any(F(fac/x) in QQ for x in allowed):
      ok=False;print('UNEXPECTED',row,odd,press,fac,flush=True)
   checks[str((row,odd,press))+'_denominator_units']=ok
lev1=load(str(OUT/'trs_B_infinity_L0_leaders_weight1.sobj'))
lev2=load(str(OUT/'trs_B_infinity_L0_leaders_wall2.sobj'))
lev3=load(str(OUT/'trs_B_infinity_L0_leaders_wall3.sobj'))
K=lev1['field'];w,p,j,P,I2,ell=K.gens();ii=K.base_ring().gen();v=w*w;Q=v-p*p
def item(lev,row,odd,press):return next(z for z in lev['results'] if z['row']==row and z['U_parity']==odd and z['pressure']==press)
def coef(lev,row,odd,press):return K(item(lev,row,odd,press)['coefficient_without_U'])
for row in [0,1]:
 entries=[z for z in lev1['results'] if z['row']==row]
 checks[str(row)+'_N_whole_regular']=all(z['order_with_U']>=row for z in entries)
 checks[str(row)+'_KC0_N_tau_divisible']=all(z['order_with_U']>=row+1 for z in entries if z['pressure'][0]==0)
 entries2=[z for z in lev2['results'] if z['row']==row]
 checks[str(row)+'_N_in_maximal_ideal']=all(z['order_with_U']+3-row>=1 for z in entries2)
 checks[str(row)+'_only_KC_KC2_in_linear_matrix']=all(z['order_with_U']+3-row>=2 for z in entries2 if z['pressure'] not in [(1,0),(2,0)])
 entries3=[z for z in lev3['results'] if z['row']==row and z['pressure'][0]==0]
 checks[str(row)+'_KC0_M_in_J']=all(z['order_with_U']+5-row>=2 for z in entries3)
 checks[str(row)+'_KB_absent_from_J_mod_mJ']=all(z['order_with_U']+5-row>=3 for z in entries3 if z['pressure'][1]>0)
 # KC-linear x-column and KC-quadratic tau-column of N, divided by KC / KC^2.
 a_even=w**2/j*p*(3*v-4*p*p)/Q
 a_odd=ii*w**3/j*(v-4*p*p)
 b_even=w**8/j**2*(v-2*p*p)/Q
 b_odd=-2*ii*p*w**9/j**2
 targets2=[(a_even,a_odd),(b_even,b_odd)]
 for press,(aa,bb) in zip([(1,0),(2,0)],targets2):
  for odd,tgt in enumerate([aa,bb]):
   power=ell if press==(1,0) else 1
   rowunit=(-ii*w if press==(1,0) else ii*w) if row else 1
   checks[str((row,odd,press))+'_KC_nonzero_matrix']=coef(lev2,row,odd,press)*ell**3==rowunit*tgt*power
 # Constant x / tau^2 columns of M=N/tau when KC=0.
 c_even=6*ii*w**3/j*(v-4*p*p)
 c_odd=6*p*w**4/j*(3*v-4*p*p)
 d_even=36*w**10/j**2*(v-2*p*p)
 d_odd=-72*ii*p*w**11*Q/j**2
 for odd,(cc,dd) in enumerate([(c_even,d_even),(c_odd,d_odd)]):
  target=cc*ell+dd if row==0 else -ii*w*cc*ell+ii*w*dd
  checks[str((row,odd))+'_KC0_J_matrix']=coef(lev3,row,odd,(0,0))*ell**3==target
checks['KC_x_unit_positive_norm']=p*p*(3*v-4*p*p)**2/Q**2+(v-4*p*p)**2/Q==v**3/Q**2
checks['KC_tau_unit_positive_norm']=(v-2*p*p)**2/Q**2+4*p*p/Q==v*v/Q**2
checks['KC0_x_unit_positive_norm']=(v-4*p*p)**2+p*p*(3*v-4*p*p)**2/Q==v**3/Q
checks['KC0_tau2_unit_positive_norm']=(v-2*p*p)**2+4*p*p*Q==v*v
checks['both_matrices_determinant_formula']=matrix(K,[[1,1],[-ii*w,ii*w]]).det()==2*ii*w
assert all(checks.values()),checks
out={'status':'PASS','check_count':len(checks),'checks':checks,'scope':'Every real pitch p, v>p^2, arbitrary real KC and KB; both T and U sheets. Exact-wall incoming arcs at every ramification and every contact order. The original normalized source ideal contains tau if KC !=0. If KC=0, the tau-divided ideal contains tau^2 and the original normalized ideal contains tau^3; division by tau is valid on every incoming Puiseux arc.'}
(OUT/'trs_B_infinity_L0_wall_reader.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
