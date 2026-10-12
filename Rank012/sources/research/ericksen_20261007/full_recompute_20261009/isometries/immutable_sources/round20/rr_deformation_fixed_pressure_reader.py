"""Exact fixed-angular-pressure condition on the complete second-order deformation cone."""
from sage.all import *
from pathlib import Path
import json,signal
OUT=Path(__file__).resolve().parent;signal.alarm(60)
d=load(str(OUT/'rr_deformation_second_order.sobj'));R=d['ring'];K=R.base_ring();E,k,o=K.gens()
r,z,lr,A0,A1,A2,B,W,H=R.gens();beta=(1-o*o-k*k)/o;gamma=(1+o*o-k*k)/o
dr=lambda f:R(f).derivative(r)+R(f).derivative(lr)/r
dz=lambda f:R(f).derivative(z)
dropA=lambda f:R({ex:cc for ex,cc in R(f).dict().items() if not any(ex[j] for j in [3,4,5])})
flux=lambda T:r*(dr(T[0,1])+dz(T[1,2])+2*T[0,1]/r)
FC=dropA(flux(d['T']));FI=dropA(flux(d['Tinv']))
BS=dropA(-d['source_residual'][0]/(2*k*r))
fixed=gamma*FC+E*k*k*beta*FI+2*k*o*r*r*BS
checks={'first_source_pressure_unchanged':dropA(flux(d['H']))==0,
 'first_inverse_pressure_unchanged':dropA(flux(-d['Cinv0']*d['H']*d['Cinv0']))==0}
groups={}
P=PolynomialRing(K,names=('B','W','H'));bb,ww,hh=P.gens()
for ex,cc in fixed.dict().items():
 key=tuple(ex[:3]);groups[key]=groups.get(key,P(0))+cc*bb**ex[6]*ww**ex[7]*hh**ex[8]
rep={'checks':checks,'fixed_pressure_obstruction_by_radial_degree':{str(key):str(f.factor()) for key,f in groups.items()},
 'scope':'Full second-order fixed-K obstruction; not actual nonlinear integrability.'}
resonant={key:P({ex:cc.subs({E:k*k}) for ex,cc in f.dict().items()}) for key,f in groups.items()}
rep['resonance_E_equals_k_squared']={str(key):str(f.factor()) if f else '0' for key,f in resonant.items()}
save({'ring':R,'FC':FC,'FI':FI,'Bstar':BS,'fixed_pressure_obstruction':fixed,'coefficient_ring':P,
 'grouped_coefficients':groups,'checks':checks},str(OUT/'rr_deformation_fixed_pressure_reader.sobj'))
(OUT/'rr_deformation_fixed_pressure_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
print(json.dumps(rep,indent=2),flush=True)
