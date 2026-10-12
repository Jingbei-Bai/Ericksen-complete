"""Exact full-base whole-source inclusion in the yw nilpotent corner ideal.

Use a literal monomial normal form, not frozen-base spreading. The
quotient relations make p31=p32=0 and every higher delta layer vanish
before original coefficient substitution. All positive powers are kept.
"""
from sage.all import *
from pathlib import Path
import time,resource,os,json
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2);p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_single_RR3_normal.sobj'));OP=d['ring'];OK=OP.base_ring();K=PolynomialRing(QQ,names=('gi_k','gi_r','gi_H','gi_Z','gi_W')).fraction_field();k,r,H,Z,W=K.gens();m=1+k;q=1+r;ph=OK.ring().hom([k/r,r,H,Z],K);T=PolynomialRing(K,names=('gi_x','gi_s','gi_t','gi_kap'));x,s,t,kap=T.gens();c=(3*k-W*W)/(W*W+k);alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));cache={};checks=[]
def ck(n,val):assert bool(val),n;checks.append(n)
def nf(f):
 out={}
 for ee,a in T(f).dict().items():
  i,j,l,n=map(int,ee)
  if j>=2 or i>=3 or (j and l+n) or (i==2 and (j or l+n)):continue
  if i==2:key=(1,1,0,0);a*=c
  elif i==1 and j==0:key=(1,0,0,l+n);a*=(-2)**l
  else:key=(i,j,l,n)
  out[key]=out.get(key,K.zero())+a
 return T(out)
def mul(a,b):return nf(a*b)
def cv(a,field,hh):
 a=field(a);key=(str(field),a)
 if key not in cache:cache[key]=hh(a.numerator())/hh(a.denominator())
 return cache[key]
def adapted(f,im):
 out=T.zero()
 for ee,a in f.dict().items():
  val=T.one()
  for v,n in zip(im,ee):
   for unused in range(int(n)):val=mul(val,v)
  if val:out+=cv(a,OK,ph)*val
 return nf(out)
graph=[adapted(f,(x,s-x,T(W),T.zero(),T.zero(),T.zero(),T.zero())) for f in d['graph_images']];es=graph[3:6];de=graph[6];rs=(alpha*W/beta*t+kap*es[1],t,T.one());ps=(x,T.one(),s-x,T(W),nf((t*(s-x)-kap*es[2])/beta),nf(x*W*t/beta+kap*(x*es[1]-es[0])/alpha));ck('literal p31 quotient zero',ps[4]==0);ck('literal p32 quotient zero',ps[5]==0);im=ps+rs
def mon(ee):
 if int(ee[4])+int(ee[5]):return T.zero()
 if int(ee[0])+int(ee[2])>=3:return T.zero()
 val=T.one()
 for v,n in zip(im,ee):
  for unused in range(int(n)):val=mul(val,v)
 return val
def norm(f,offset,discard,field,hh):
 ans=T.zero()
 for ee,a in f.dict().items():
  pd,rd=sum(ee[:6]),sum(ee[6:])
  if (pd,rd) in discard:continue
  ke,dp=offset[0]-2*pd-3*rd,offset[1]-pd-rd;ck('full layer nonnegative '+str((offset,pd,rd)),ke>=0 and dp>=0)
  if dp>=2 or (dp and ke):continue
  val=mon(ee)
  if dp:val=mul(val,de)
  if ke:val=mul(val,kap**ke)
  if val:ans+=cv(a,field,hh)*val
 return nf(ans)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));R=qf['ring'];sources=[];DC=matrix(R,qf['D'])*matrix(R,qf['C'])
for j,full in enumerate(qf['full_R_sources']):
 ans=norm(full,(8,4),{(3,1)},OK,ph)+sum((mul(adapted(DC[j,i],im),es[i]) for i in range(3)),T.zero());ck('complete R in generic ideal '+str(j),nf(ans)==0);sources.append(nf(ans))
ans=norm(qf['full_Q0'],(7,3),{(1,2)},OK,ph)+sum((mul(adapted(qf['Q0_L_cofactors'][i],im),es[i]) for i in range(3)),T.zero());ck('complete Q0 in generic ideal',nf(ans)==0);sources.append(nf(ans));print('R Q PASS',time.monotonic()-st,flush=True)
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'));old=next(iter(hs['items'][0]['H41_original'].values())).parent();oph=old.ring().hom([k/r,k*q/(r*m),H/r,Z*q/r],K)
def orig(dd):
 ans=T.zero()
 for ee,a in dd.items():
  val=mon(ee)
  if val:ans+=cv(a,old,oph)*val
 return nf(ans)
for j,it in enumerate(hs['items']):
 ans=sum((mul(orig(it['L_cofactors_original'][i]),es[i]) for i in range(3)),T.zero())
 # Complete weight bookkeeping leaves only H41, kappa*H50,
 # delta*H13 besides the replaced H32 module. Every other layer
 # contains delta^2 or kappa*delta, both literally zero here.
 ans+=orig(it['H41_original'])
 h13=d['boundary'][7+j].derivative(OP.gen(6))
 ans+=mul(de,adapted(h13,(x,s-x,T(W),T.zero(),T.zero(),T.zero(),T.zero())))
 h50=d['kappa_native_extra'][7+j]
 ans+=mul(kap,adapted(h50,(x,s-x,T(W),T.zero(),T.zero(),T.zero(),T.zero())))
 ck('complete H in generic ideal '+str(j),nf(ans)==0);sources.append(nf(ans));print('H',j,'PASS',time.monotonic()-st,flush=True)
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));ell=matrix(R,qf['L'])*vector(R,R.gens()[6:])
for it in ev['items']:
 rem=it['full_source']-OK(it['square_unit'])*ell[int(it['square_component'])]**2-sum((R(it['A_p_quadratic'][i])*ell[i] for i in range(3)),R.zero());ans=norm(rem,(8,4),set(),OK,ph)+cv(it['square_unit'],OK,ph)*mul(es[int(it['square_component'])],es[int(it['square_component'])])+sum((mul(adapted(it['A_p_quadratic'][i],im),es[i]) for i in range(3)),T.zero());ck('complete even in generic ideal '+str(it['Q_index']),nf(ans)==0);sources.append(nf(ans))
save({'ring':T,'graph':graph,'chart_p':ps,'chart_Rbar':rs,'ideal_generators':(kap*s,t*s,s*s,(t+2*kap)*x,x*x-c*x*s),'source_normal_forms':sources,'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_generic_relative_ideal.sobj'),compress=False);(p/'independent_subcritical_single_RR_yw_generic_relative_ideal.json').write_text(json.dumps({'status':'PASS','sources':len(sources),'checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS ALL28 generic relative ideal',time.monotonic()-st,flush=True)
