"""Exact whole cyclic-base 31source inclusion, with k and w independent at the minus attachment."""
from sage.all import *
from pathlib import Path
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_single_RR3_normal.sobj'));pa=load(str(p/'independent_subcritical_single_RR_full_pair_chart.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2);print('PID',os.getpid(),flush=True);OP=d['ring'];OK=OP.base_ring();K=PolynomialRing(QQ,names=('mi_k','mi_r','mi_W')).fraction_field();k,r,W=K.gens();H,Z=K(1),K(-1);m=1+k;q=1+r;ph=OK.ring().hom([k/r,r,H,Z],K);T=PolynomialRing(K,names=('mi_x','mi_s','mi_t','mi_kap'));x,s,t,kap=T.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));AA=2*m*(k+W*W)/k;cache={};checks=[]
def ck(n,v):assert bool(v),n;checks.append(n)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def nf(f):
 out={}
 for ee,a in T(f).dict().items():
  i,j,l,n=map(int,ee)
  if j>=2 or (j and l+n) or (i>=3 and l+n):continue
  if i>=1 and l:a*=(-2)**l;n+=l;l=0
  key=(i,j,l,n);out[key]=out.get(key,K.zero())+a
 clean={}
 for ee,a in out.items():
  a=nk(a)
  if ee[0]==2 and ee[3]>=1 and a:
   DD=K.ring()(k+W*W);rem=a.numerator().quo_rem(DD)[1]
   if not rem and a.denominator().quo_rem(DD)[1]:continue
  if a:clean[ee]=a
 return T(clean)
def mul(a,b):return nf(a*b)
def cv(a,field,hh):
 a=field(a);key=(str(field),a)
 if key not in cache:cache[key]=nk(hh(a.numerator())/hh(a.denominator()))
 return cache[key]
def adapted(f,im,field=OK,hh=ph):
 out=T.zero()
 for ee,a in f.dict().items():
  val=T.one()
  for v,n in zip(im,ee):
   for unused in range(int(n)):val=mul(val,v)
   if not val:break
  if val:out+=cv(a,field,hh)*val
 return nf(out)
graph=[adapted(f,(x,s-x,T(W),T.zero(),T.zero(),T.zero(),T.zero())) for f in d['graph_images']];es=graph[3:6];de=graph[6];rs=(alpha*W/beta*t+kap*es[1],t,T.one());ps=(x,T.one(),s-x,T(W),nf((t*(s-x)-kap*es[2])/beta),nf(x*W*t/beta+kap*(x*es[1]-es[0])/alpha));ck('literal p31 zero',ps[4]==0);ck('literal p32 zero',ps[5]==0);im=ps+rs
def mon(ee):
 if int(ee[4])+int(ee[5]):return T.zero()
 val=T.one()
 for v,n in zip(im,ee):
  for unused in range(int(n)):val=mul(val,v)
  if not val:break
 return val
def norm(f,offset,discard,field,hh):
 ans=T.zero()
 for ee,a in f.dict().items():
  pd,rd=sum(ee[:6]),sum(ee[6:])
  if (pd,rd) in discard:continue
  ke,dp=offset[0]-2*pd-3*rd,offset[1]-pd-rd;ck('layer regular '+str((offset,pd,rd)),ke>=0 and dp>=0)
  if dp>=2 or (dp and ke):continue
  val=mon(ee)
  if dp:val=mul(val,de)
  if ke:val=mul(val,kap**ke)
  if val:ans+=cv(a,field,hh)*val
 return nf(ans)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));R=qf['ring'];sources=[];DC=matrix(R,qf['D'])*matrix(R,qf['C'])
for j,f in enumerate(qf['full_R_sources']):
 ans=norm(f,(8,4),{(3,1)},OK,ph)+sum((mul(adapted(DC[j,i],im),es[i]) for i in range(3)),T.zero());ck('full R '+str(j),nf(ans)==0);sources.append(nf(ans))
ans=norm(qf['full_Q0'],(7,3),{(1,2)},OK,ph)+sum((mul(adapted(qf['Q0_L_cofactors'][i],im),es[i]) for i in range(3)),T.zero());ck('full Q0',nf(ans)==0);sources.append(nf(ans));print('RQ PASS',time.monotonic()-st,flush=True)
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'));old=next(iter(hs['items'][0]['H41_original'].values())).parent();oph=old.ring().hom([k/r,k*q/(r*m),H/r,Z*q/r],K)
def orig(dd):
 ans=T.zero()
 for ee,a in dd.items():
  val=mon(ee)
  if val:ans+=cv(a,old,oph)*val
 return nf(ans)
for j,it in enumerate(hs['items']):
 ans=orig(it['H41_original'])+sum((mul(orig(it['L_cofactors_original'][i]),es[i]) for i in range(3)),T.zero());ans+=mul(de,adapted(d['boundary'][7+j].derivative(OP.gen(6)),(x,s-x,T(W),T.zero(),T.zero(),T.zero(),T.zero())));ans+=mul(kap,adapted(d['kappa_native_extra'][7+j],(x,s-x,T(W),T.zero(),T.zero(),T.zero(),T.zero())));ck('full H '+str(j),nf(ans)==0);sources.append(nf(ans));print('H',j,'PASS',time.monotonic()-st,flush=True)
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));ell=matrix(R,qf['L'])*vector(R,R.gens()[6:])
for it in ev['items']:
 rem=it['full_source']-OK(it['square_unit'])*ell[int(it['square_component'])]**2-sum((R(it['A_p_quadratic'][i])*ell[i] for i in range(3)),R.zero());ans=norm(rem,(8,4),set(),OK,ph)+cv(it['square_unit'],OK,ph)*mul(es[int(it['square_component'])],es[int(it['square_component'])])+sum((mul(adapted(it['A_p_quadratic'][i],im),es[i]) for i in range(3)),T.zero());ans=nf(ans)
 if int(it['Q_index'])==6:
  C=adapted(d['graph_boundary'][27],(x,s-x,T(W),T.zero(),T.zero(),T.zero(),T.zero()));ck('full even6 remainder exactly boundary C',ans==C);print('EXACT C',C,'FACT',C.factor(),flush=True)
 else:ck('full even '+str(it['Q_index']),ans==0)
 sources.append(ans)
PK=pa['ring'].base_ring();pph=PK.ring().hom([k,r,H,Z],K);pims=(x,s-x,T(W))+tuple(es)+(de,t,kap)
for it in pa['items']:
 ans=adapted(it['source'],pims,PK,pph);ck('full pair '+str(it['pair']),ans==0);sources.append(ans)
save({'ring':T,'graph':graph,'chart_p':ps,'chart_Rbar':rs,'ideal_generators':(kap*s,t*s,s*s,(t+2*kap)*x,C,kap*x**3),'source_normal_forms':sources,'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_minus_cyclic_relative_ten_ideal.sobj'),compress=False);(p/'independent_subcritical_single_RR_yw_minus_cyclic_relative_ten_ideal.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'sources':len(sources),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS ALL31 CYCLIC TEN',time.monotonic()-st,flush=True)
