"""Complete original sources, first normal along whole horizontal family.
All H,Z,t,v,kappa are retained. Only the explicitly named transverse
normal ideal is squared out; this is an exact conormal computation.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic()
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));R=qf['ring'];OK=R.base_ring()
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'));raw=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));old=load(str(p/'independent_delta_H_global_L_module.sobj'))['original_field'];ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
wb=load(str(p/'pure_x_root_generic_w_bridge.sobj'));wf=load(str(p/'pure_x_root_generic_w_reader.partial.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
print('PID',os.getpid(),flush=True);K=PolynomialRing(QQ,names=('k','r','H','Z')).fraction_field();k,r,H,Z=K.gens();m=1+k;q=1+r;P=PolynomialRing(K,names=('y','z','w','e0','e1','e2','delta','t','v','kap'));y,z,w,e0,e1,e2,de,t,v,kap=P.gens();normal=tuple(range(7));ph=OK.ring().hom([k/r,r,H,Z],K);oph=old.ring().hom([k/r,k*q/(r*m),H/r,q*Z/r],K);cache={};ocache={};checks=[]
def ck(n,x):assert bool(x),n;checks.append(n)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=nk(ph(c.numerator())/ph(c.denominator()))
 return cache[c]
def oc(c):
 c=old(c)
 if c not in ocache:ocache[c]=nk(oph(c.numerator())/oph(c.denominator()))
 return ocache[c]
def nf(f):return P({tuple(map(int,e)):nk(c) for e,c in f.dict().items() if sum(e[:7])<=1 and c})
al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));shift=-2*H*(Z+1)/(Z*(H+1));oldv=(al+kap*shift)*v;olde0=e0+shift*v
# Here e0 names the true transverse E, and v names the horizontal l.
imgs=(P.one(),y,z,w,(t*z-kap*e2)/be,v-kap*e0/al,oldv,t,P.one());es=(olde0,e1,e2);inc=oldv*y-al/be*t*w-kap*e1
power={}
def mon(ee):
 ans=P.one()
 for j,n in enumerate(ee):
  if not n:continue
  ky=(j,int(n))
  if ky not in power:power[ky]=nf(imgs[j]**int(n))
  ans=nf(ans*power[ky])
  if not ans:break
 return ans
def mp(f,mapper=cv):
 ans=P.zero()
 for ee,c in f.dict().items():
  mm=mon(ee)
  if mm:ans=nf(ans+mapper(c)*mm)
 return ans
def norm(f,offset,discard):
 ans=P.zero()
 for ee,c in f.dict().items():
  pd=sum(ee[:6]);rd=sum(ee[6:]);ee=tuple(map(int,ee))
  if (pd,rd) in discard:continue
  ke,dd=offset[0]-2*pd-3*rd,offset[1]-pd-rd;assert ke>=0 and dd>=0
  if dd>1:continue
  mm=nf(de**dd*mon(ee))
  if mm:ans=nf(ans+cv(c)*kap**ke*mm)
 return ans
DC=matrix(R,qf['D'])*matrix(R,qf['C']);sources=[]
for j,rf in enumerate(qf['full_R_sources']):sources.append(nf(norm(rf,(8,4),{(3,1)})+sum((nf(mp(DC[j,i])*es[i]) for i in range(3)),P.zero())))
sources.append(nf(norm(qf['full_Q0'],(7,3),{(1,2)})+sum((nf(mp(qf['Q0_L_cofactors'][i])*es[i]) for i in range(3)),P.zero())))
print('R/Q',time.monotonic()-st,flush=True);pool=raw['polynomial_pool']
for j,it in enumerate(hs['items']):
 acc={}
 for sign,lr,rr in raw['records'][j]['full_signed_product_DAG']:
  for le,lc in pool[lr].items():
   for re,rc in pool[rr].items():
    ee=tuple(int(a)+int(b) for a,b in zip(le,re));pd=sum(ee[:6]);rd=sum(ee[6:])
    if (pd,rd)==(3,2):continue
    ke,dd=11-2*pd-3*rd,5-pd-rd
    if dd>1:continue
    # Original degree-six terms can occur in the signed DAG and cancel.
    if ke>=0 and dd>=0 and not nf(de**dd*mon(ee)):continue
    acc[ee]=nk(acc.get(ee,K.zero())+sign*oc(lc)*oc(rc))
 ans=P.zero()
 for i in range(3):
  dd=it['L_cofactors_original'][i];tmp=P.zero()
  for ee,c in dd.items():
   mm=nf(es[i]*mon(ee))
   if mm:tmp=nf(tmp+oc(c)*mm)
  ans=nf(ans+tmp)
 for ee,c in acc.items():
  if not c:continue
  pd=sum(ee[:6]);rd=sum(ee[6:]);ke,dd=11-2*pd-3*rd,5-pd-rd;ck('complete DAG layer bounds '+str((j,pd,rd)),ke>=0 and dd>=0)
  ans=nf(ans+c*kap**ke*de**dd*mon(ee))
 sources.append(ans);print('H',j,'normal terms',len(ans.dict()),flush=True)
def cp(f):return R({tuple(map(int,e))+(0,)*(9-len(e)):OK(c) for e,c in f.dict().items()})
ell=matrix(R,qf['L'])*vector(R,R.gens()[6:])
for it in ev['items']:
 src=cp(it['full_source']);rem=src-OK(it['square_unit'])*ell[int(it['square_component'])]**2-sum((cp(it['A_p_quadratic'][i])*ell[i] for i in range(3)),R.zero())
 ans=nf(norm(rem,(8,4),set())+cv(it['square_unit'])*es[int(it['square_component'])]**2+sum((nf(mp(cp(it['A_p_quadratic'][i]))*es[i]) for i in range(3)),P.zero()));sources.append(ans)
shift_images=(y,z,w,olde0,e1,e2,de,t,oldv,kap);shift_powers={}
def shiftpull(f):
 ans=P.zero()
 for ee,c in f.dict().items():
  val=P(c)
  for j,n in enumerate(ee):
   if not n:continue
   ky=(j,int(n))
   if ky not in shift_powers:shift_powers[ky]=nf(shift_images[j]**int(n))
   val=nf(val*shift_powers[ky])
   if not val:break
  ans=nf(ans+val)
 return ans
gs=[shiftpull(P(f)) for f in wf['full_G_sources']];vec=wb['combination'];wn=nf(vec[0]*gs[0]+(vec[1]*w+vec[2]*oldv)*gs[1]+vec[3]*y*gs[2]-shiftpull(P(wb['incidence_cofactor']))*inc)
ck('whole normal W numerator kappa divisible',all(e[9]>=1 for e in wn.dict()));W=P({tuple(int(a)-(1 if i==9 else 0) for i,a in enumerate(e)):c for e,c in wn.dict().items()});sources.append(W)
for i,f in enumerate(sources):ck('whole source horizontal zero '+str(i),all(sum(e[:7])==1 for e in f.dict()))
save({'ring':P,'sources':sources,'W_index':28,'normal_columns':normal,'shift_images':shift_images,'shift':shift,'checks':checks,'scope':__doc__},str(p/'pure_x_root_generic_horizontal_normals.sobj'),compress=False)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__};(p/'pure_x_root_generic_horizontal_normals.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
