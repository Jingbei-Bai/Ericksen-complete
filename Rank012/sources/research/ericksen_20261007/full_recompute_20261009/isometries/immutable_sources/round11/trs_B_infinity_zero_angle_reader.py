"""Assemble all-pitch whole L->0 coverage from original-source coefficients."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent
lev=load(str(OUT/'trs_B_infinity_L0_leaders_weight1.sobj'))
new=load(str(OUT/'trs_B_infinity_L0_Newton_reader.sobj'))
K=lev['field'];w,p,j,P,I2,ell=K.gens();ii=K.base_ring().gen();v=w*w;Q=v-p*p;D=w-ii*ell
checks={}
def term(odd,press):return next(z for z in new['leaders'] if z['U_parity']==odd and z['pressure']==press)
for odd in [0,1]:
 for press,tp,lp in [((0,0),2,0),((0,1),3,0),((1,0),2,-1),((1,1),3,-1),((2,0),4,-2)]:
  z=term(odd,press);checks[str((odd,press))+'_Newton_weight']=z['tau_power_with_U']==tp and z['L_power']==lp
z=term(0,(0,2));checks['KB2_Newton_weight']=z['tau_power_with_U']==5 and z['L_power']==-1
I1=(v**3+v*I2-j-P)/v**2
for odd,target in enumerate([-14*p*v/j,-7*ii*w**3*(v-2*p*p)/j]):
 checks[str(odd)+'_Newton_KC_coefficient']=K(term(odd,(1,0))['coefficient_without_U_unit'])==target
for odd,target in enumerate([6*I1*(v-2*p*p)/j,-12*ii*I1*p*w*Q/j]):
 checks[str(odd)+'_Newton_constant_coefficient']=K(term(odd,(0,0))['coefficient_without_U_unit'])==target
checks['Newton_KC_coefficient_norm']=4*p*p+(v-2*p*p)**2/Q==v*v/Q
checks['Newton_constant_coefficient_norm']=(v-2*p*p)**2+4*p*p*Q==v*v
def walloff(odd,press):return next(z for z in lev['results'] if z['row']==0 and z['U_parity']==odd and z['pressure']==press)
targets={(0,(1,0)):-p*w*w*(3*v-4*p*p)/(j*Q*D*D),
 (1,(1,0)):-ii*w**3*(v-4*p*p)/(j*D*D),
 (0,(0,0)):-6*ii*w**3*(v-4*p*p)/(j*D*D),
 (1,(0,0)):-6*p*w**4*(3*v-4*p*p)/(j*D*D)}
for (odd,press),tgt in targets.items():
 z=walloff(odd,press);checks[str((odd,press))+'_offwall_exact_coefficient']=K(z['coefficient_without_U'])==tgt
checks['offwall_KC_unique_lowest']=all(z['order_with_U']>=1 for z in lev['results'] if z['row']==0 and z['pressure']!=(1,0))
checks['offwall_KC0_constant_unique_lowest']=all(z['order_with_U']>=2 for z in lev['results'] if z['row']==0 and z['pressure'][0]==0 and z['pressure'][1]>0)
checks['offwall_KC_norm']=p*p*(3*v-4*p*p)**2/Q**2+(v-4*p*p)**2/Q==v**3/Q**2
checks['offwall_constant_norm']=(v-4*p*p)**2+p*p*(3*v-4*p*p)**2/Q==v**3/Q
assert all(checks.values()),checks
parts={name:json.loads((OUT/name).read_text()) for name in ['trs_B_infinity_L0_Newton_reader.json','trs_B_infinity_L0_wall_reader.json']}
assert all(x['status']=='PASS' for x in parts.values())
out={'status':'PASS','check_count':len(checks)+sum(x['check_count'] for x in parts.values()),'assembly_checks':checks,'source_readers':{name:x['check_count'] for name,x in parts.items()},'scope':'All B-infinity arcs with L->0, every real p, v>p^2, positive material spectrum, every real KC and KB; every Newton slope and every exact-wall contact order. Both quadratic sheets by complex conjugation and U-sign.'}
(OUT/'trs_B_infinity_zero_angle_reader.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
