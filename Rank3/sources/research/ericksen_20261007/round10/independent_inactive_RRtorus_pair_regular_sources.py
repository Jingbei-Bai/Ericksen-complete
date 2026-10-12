from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));lf=load(str(p/'independent_inactive_RRtorus_lifted_fields_cyclic.sobj'));Q=lf['ring'];K=Q.base_ring();a,r=K.gens();m=1+a*r;q=1+r;F=Q.fraction_field();l,rho,k,s,w,e0,e1,e2,de=Q.gens();K0=d['ring'].base_ring();ph=K0.ring().hom((a,a*q/m,1/r,-q/r),K);cache={}
def cf(c):
 c=K0(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
pv=[F(1),(s*rho+k*e1)/(l+k*e0),(w+k*e2)/rho,s,w,l];Rv=[l+k*e0,rho,F(1)];out=[];ns=[s,w,e0,e1,e2,de];sub={v:0 for v in ns}
for it in d['items']:
 layers={}
 for ex,co in it['full'].dict().items():
  ex=tuple(map(int,ex));ii=sum(ex[:6]);jj=sum(ex[6:]);ke=8-2*ii-3*jj;dd=4-ii-jj
  assert dd>=0
  z=F(cf(co))*prod(F(pv[j])**ex[j] for j in range(6))*prod(Rv[j]**ex[j+6] for j in range(3))*de**dd
  layers[ke]=layers.get(ke,F(0))+z
 ff=sum((z*k**ke for ke,z in layers.items()),F(0));pole=0
 while ff.denominator().subs({k:0})==0:
  ff*=k;pole+=1
 assert pole<=2
 print('POLE_REMOVED',pole,flush=True)
 assert ff.subs(sub)==0
 print('SOURCE',it['pair'],'TERMS',len(ff.numerator().dict()),flush=True)
 fk=ff.subs({k:0});N=vector(F,[fk.derivative(v).subs(sub) for v in ns]);out.append(dict(pair=it['pair'],kappa_power=pole,full_source=ff,kappa_zero=fk,normal=N))
 save(dict(ring=Q,items=out,scope='Full three global pair sources, weight6 homogenization divided by tau squared and multiplied by the recorded minimal kappa power, cyclic all materials; regularity and whole S-prime zero checked literally.'),str(p/'independent_inactive_RRtorus_pair_regular_sources.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
