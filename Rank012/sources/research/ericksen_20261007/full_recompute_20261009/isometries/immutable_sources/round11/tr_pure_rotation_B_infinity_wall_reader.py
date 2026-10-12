"""Whole-ring bridge for the entire near-wall B=infinity incoming stratum.

This reader checks the required denominator support and source-linked exact
jets. Analyticity makes membership in (x,tau^2) a FINITE quotient check;
no claim that an arbitrary leading ideal determines a whole source is made.
"""
from sage.all import *
from pathlib import Path
import json,resource
OUT=Path(__file__).resolve().parent
dat=load(str(OUT/'tr_pure_rotation_first_order_fields.sobj'));P=dat['P'];K=P.fraction_field()
v,j,A,T,B,L=map(K,P.gens())
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
checks={}
sources=[load(str(OUT/('tr_pure_rotation_'+lab+'_frobenius.sobj')))['coefficients'] for lab in ['logU','L']]
factorlist=list(sources[1]['1'].denominator().factor())
H=max((ff for ff,ee in factorlist),key=lambda ff:len(ff.dict()))
checks['wall_polynomial_linear_L']=H.degree(P.gen(5))==1
for row,cols in enumerate(sources):
 for label,f in cols.items():
  rem=K(f.denominator())/K(H)**3
  ok=True
  for poly in [rem.numerator(),rem.denominator()]:
   for ff,ee in poly.factor():
    if ff not in [P.gen(0),P.gen(1),P.gen(2),P.gen(3),P.gen(4),P.gen(0)+P.gen(5)**2]:ok=False
  checks['row'+str(row)+'_'+label+'_only_wall_cubed_and_units']=ok
lev1=load(str(OUT/'tr_pure_rotation_B_infinity_leaders_Lweight1.sobj'))['leaders']
lev2=load(str(OUT/'tr_pure_rotation_B_infinity_wall_shift_2.sobj'))['leaders']
lev3=load(str(OUT/'tr_pure_rotation_B_infinity_wall_shift_3.sobj'))['leaders']
for row,name in enumerate(['logU_frobenius','L_frobenius']):
 shift=0 if row==0 else 1
 checks[name+'_N_regular_all_columns']=all(int(z['order_with_E'])>=shift for z in lev1[name].values())
 checks[name+'_KC0_N_tau_divisible']=all(int(lev1[name][lab]['order_with_E'])>=1+shift for lab in ['1','ED','DD'])
 # x=ell*tau^2. M1=x^3/tau F1; M2=x^3/tau^2 F2.
 checks[name+'_KC0_M_in_x_tau2']=all(int(lev3[name][lab]['order_with_E'])+5-row>=2 for lab in ['1','ED','DD'])
 coeff=lev3[name]['1']['coefficient_without_E']; KK=coeff.parent();ww,jj,pp,ii2,ell=KK.gens();II=KK.base_ring().gen()
 target=(6*II*ww**5/jj*ell+36*ww**12/jj**2) if row==0 else (6*ww**6/jj*ell+36*II*ww**13/jj**2)
 checks[name+'_KC0_J2_row']=coeff*ell**3==target
 # At x=ell*tau, N_i have total degree-one coefficients KC*x,KC^2*tau.
 ec=lev2[name]['EC']['coefficient_without_E'];cc=lev2[name]['CC']['coefficient_without_E'];K2=ec.parent();ww2,jj2,pp2,ii22,ell2=K2.gens();II2=K2.base_ring().gen()
 targetec=II2*ww2**3/jj2*ell2 if row==0 else ww2**4/jj2*ell2
 targetcc=ww2**8/jj2**2 if row==0 else II2*ww2**9/jj2**2
 checks[name+'_KC_nonzero_x_row']=ec*ww2**2*ell2**3==targetec
 checks[name+'_KC_nonzero_tau_row']=cc*ell2**3==targetcc
K3=lev3['logU_frobenius']['1']['coefficient_without_E'].parent();ww,jj,pp,ii2,ell=K3.gens();II=K3.base_ring().gen()
det0=(6*II*ww**5/jj)*(36*II*ww**13/jj**2)-(6*ww**6/jj)*(36*ww**12/jj**2)
detC=(II*ww**3/jj)*(II*ww**9/jj**2)-(ww**4/jj)*(ww**8/jj**2)
checks['KC0_J2_unit_determinant']=det0==-432*ww**18/jj**3
checks['KC_nonzero_linear_unit_determinant_without_KC3']=detC==-2*ww**12/jj**3
assert all(checks.values()),checks
out={'status':'PASS','check_count':len(checks),'checks':checks,'scope':'All B-infinity arcs with L/tau tending to the exact finite determinant-wall leading root; arbitrary ramification and both root signs. KC=0/nonzero, every KB retained. No flatness source used.'}
(OUT/'tr_pure_rotation_B_infinity_wall_reader.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
