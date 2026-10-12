"""Exact full Q6 own-active G regularity on the RR3,p13 charts.

Build the two genuinely dangerous original-source numerator layers, with
every base coefficient derivative. Verify L² and L membership in the
exact smooth incidence chart with free e and t, before conic restriction.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();g=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));it=g['items'][1];R=ev['ring'];K=R.base_ring();a,r,H,Z=K.gens();q=1+r;m=1+a*r;v=R.gens();old=g['ring'].base_ring();ph=old.ring().hom([a,a*q/m,H/r,Z*q/r],K);cache={};checks=[]
def ck(label,value):
 assert bool(value),label;checks.append(label)
def cv(c):
 c=old(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def cp(f):return R({tuple(map(int,ex)):cv(c) for ex,c in f.dict().items()})
fields=[cp(f) for f in it['field_numerators']];bases=[cp(f) for f in it['base_numerators']];Q=next(ee['full_source'] for ee in ev['items'] if int(ee['Q_index'])==6);parts=[(fields[j],Q.derivative(v[j])) for j in range(9)]
for j,param in enumerate((H,Z)):
 scale=r if j==0 else r/q;der=R({tuple(map(int,ex)):scale*c.derivative(param) for ex,c in Q.dict().items() if c.derivative(param)});parts.append((bases[j],der))
def face(f,pd,rd):return R({tuple(map(int,ex)):c for ex,c in f.dict().items() if sum(ex[:6])==pd and sum(ex[6:])==rd})
layers={rd:R.zero() for rd in (2,3)}
for left,right in parts:
 ll=sorted(set((sum(ex[:6]),sum(ex[6:])) for ex in left.dict()));rr=sorted(set((sum(ex[:6]),sum(ex[6:])) for ex in right.dict()))
 for li,lj in ll:
  for ri,rj in rr:
   if li+lj+ri+rj==6 and lj+rj in layers:layers[lj+rj]+=face(left,li,lj)*face(right,ri,rj)
print('LAYERS',[len(f.dict()) for f in layers.values()],'sec',time.monotonic()-st,flush=True)
P=PolynomialRing(K,names=('rg_x','rg_z','rg_w','rg_t','rg_e0','rg_e1','rg_e2','rg_kappa'));x,z,w,t,e0,e1,e2,kap=P.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));images=(x,P.one(),z,w,(t*z-kap*e2)/beta,x*w*t/beta+kap*(x*e1-e0)/alpha,alpha/beta*w*t+kap*e1,t,P.one());hm=R.hom(images,P)
ells=ev['L'].change_ring(R)*vector(R,v[6:]);ck('ell0 exact chart',hm(ells[0])==kap*e0);ck('ell1 exact chart',hm(ells[1])==kap*e1);ck('ell2 exact chart',hm(ells[2])==kap*e2)
mapped={};quotients={}
for rd,power in ((3,2),(2,1)):
 mapped[rd]=hm(layers[rd]);ck('complete dangerous G6 RR%d lies in ell^%d'%(rd,power),all(ex[7]>=power for ex in mapped[rd].dict()));quotients[rd]=P({tuple(int(n)-(power if j==7 else 0) for j,n in enumerate(ex)):c for ex,c in mapped[rd].dict().items()});ck('literal exact kappa division G6 RR'+str(rd),mapped[rd]==kap**power*quotients[rd]);print('CHECK',rd,'quotient terms',len(quotients[rd].dict()),'sec',time.monotonic()-st,flush=True)
ck('gauge denominator p13 is chart unit',hm(cp(it['denominator']))==1)
payload={'original_ring':R,'ring':P,'layers':layers,'mapped_layers':mapped,'exact_kappa_quotients':quotients,'chart_images':images,'gauge_item':1,'Q_index':6,'checks':checks,'scope':__doc__,'normalization':'G_h/(tau² sigma), tau=kappa²delta, sigma=kappadelta','full_Q_provenance':'independent_delta_subcritical_even_full_source.sobj Q6 exact source'}
save(payload,str(p/'independent_subcritical_single_RR_complete_E_regularity.sobj'),compress=False);(p/'independent_subcritical_single_RR_complete_E_regularity.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('DONE',len(checks),time.monotonic()-st,flush=True)
