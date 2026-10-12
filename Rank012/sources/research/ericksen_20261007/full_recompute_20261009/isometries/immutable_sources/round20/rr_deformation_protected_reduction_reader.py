"""All-parameter second-order reduction without exceptional-state echelon pivots."""
from sage.all import *
from pathlib import Path
import json,signal
OUT=Path(__file__).resolve().parent;signal.alarm(50)
d=load(str(OUT/'rr_deformation_second_order.sobj'));R=d['ring'];K=R.base_ring();E,k,w=K.gens()
r,z,lr,A0,A1,A2,B,W,H=R.gens();obs=d['obstructions'];A=A0+A1*z+A2*z*z
Delta=E*E-E*(1+w*w+k*k)+k*k
disc=(k*k+w*w-1)**2+4*w*w
def by_radial(f,n):
 return sum(R({ex:cc}) for ex,cc in R(f).dict().items() if ex[0]==n)
checks={'leading_radius_obstruction':by_radial(obs[1],-7)==-12*disc*A*A/(w**3*r**7)}
beta=(1-w*w-k*k)/w
checks['hand_plane_theta_second_term']=by_radial(d['T'][0,1],-4)==-k*(4+beta*beta)*A*A/(2*w*r**4)
checks['hand_source_angular_lowest_term']=by_radial(d['source_residual'][0],-5)==-4*k*(4+beta*beta)*A*A/(w*r**5)
small=[R({ex:cc for ex,cc in f.dict().items() if not any(ex[j] for j in [3,4,5])}) for f in obs]
checks['logarithmic_shear_protected_obstruction']=small[5]==8*E*(k*k-E)*H*H/(Delta*r**5)
Z=3*(E-1)*(k*k-1)-(E+1)*w*w
T=(E-1)*(k*k-E)*(k*k-1)+w*w*((5*E+1)*k*k-E*(E+5))
P=(E-1)*Z*B*B+4*E*w*Z*B*W+4*E*T*W*W
checks['remaining_conic_protected_source']=by_radial(small[0],1)==-k*P*r/(Delta*w)
# Verify every complete reduced coefficient is a combination of P and (E-k²)H².
multipliers=[]
for i,f in enumerate(small):
 for ex,cc in f.dict().items():
  pass
 # Distinct pressure/shear quadratic supports make a canonical projection possible.
 Q=PolynomialRing(K,names=('B','W','H'));bb,ww,hh=Q.gens()
 groups={}
 for ex,cc in f.dict().items():
  key=tuple(ex[:3]);groups[key]=groups.get(key,Q(0))+cc*bb**ex[6]*ww**ex[7]*hh**ex[8]
 p=Q((E-1)*Z)*bb*bb+4*E*w*Z*bb*ww+4*E*T*ww*ww
 for key,g in groups.items():
  cp=g.monomial_coefficient(bb**2)/((E-1)*Z) if g.monomial_coefficient(bb**2) else K(0)
  rest=g-cp*p;ch=rest.monomial_coefficient(hh**2)/(E-k*k)
  checks['reduced_source_%d_%s'%(i,key)]=rest==ch*(E-k*k)*hh**2
  # Crucial: the coefficients must not acquire E=1, Z=0, E=k² denominators.
  for nam,vv in [('P',cp),('H',ch)]:
   den=vv.denominator();facts=list(den.factor())
   protected=all(ff in [E,E.parent()(k),E.parent()(w),Delta] or any(ff==c*g0 for c in [1,-1] for g0 in [E,k,w,Delta]) for ff,power in facts)
   # Factor normalization can be nonmonic; use exact divisibility by the guard.
   guard=E*k*w*Delta
   protected=all(K(guard).numerator()%ff==0 for ff,power in facts)
   checks['protected_multiplier_%d_%s_%s'%(i,key,nam)]=bool(protected)
  multipliers.append((i,key,cp,ch))
checks['conic_discriminant_identity']=(4*E*w*Z)**2-4*(E-1)*Z*(4*E*T)==16*E*Delta*( (E-1)*(k*k-1)+(E+1)*w*w )*Z
assert all(checks.values()),{k:v for k,v in checks.items() if not v}
rep={'status':'FULL_NONZERO_TWIST_SECOND_ORDER_OBSTRUCTIONS_PROTECTED',
 'checks':checks,'check_count':len(checks),'conic':str(P),
 'guards':'E>0,k!=0,omega!=0,Delta!=0; no division by E-1, E-k² or Z',
 'conclusion':'A0=A1=A2=0; (E-k²)Wlog²=0; P(B1,W2)=0. This is formal order-two integrability, not an actual solution or full rigidity.'}
(OUT/'rr_deformation_protected_reduction_reader.json').write_text(json.dumps(rep,indent=2),encoding='utf8')
save({'checks':checks,'P':P,'Delta':Delta,'positive_discriminant':disc,'multipliers':multipliers},str(OUT/'rr_deformation_protected_reduction_reader.sobj'))
print(json.dumps(rep,indent=2),flush=True)
