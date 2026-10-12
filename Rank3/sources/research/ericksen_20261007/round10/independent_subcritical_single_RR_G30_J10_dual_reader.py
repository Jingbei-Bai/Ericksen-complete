"""Whole-base complete G30 dual tests needed by the J10 saturation proof."""
from sage.all import *
from pathlib import Path
import resource,os,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_G30_original_generic.sobj'));b=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_exact.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2);print('PID',os.getpid(),flush=True);O=d['ring'];OK=O.base_ring();K=PolynomialRing(QQ,names=('gd_k','gd_r','gd_W','gd_H','gd_Z')).fraction_field();k,r,W,H,Z=K.gens();m=1+k;q=1+r;ph=OK.ring().hom([k,r,H,Z],K);T=PolynomialRing(K,names=('gd_x','gd_s','gd_t','gd_kap'));x,s,t,kap=T.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));es=(-2*x*W+(4*H*Z-2*H-4*Z**2+2*Z)/(H*Z)*s*W,-2*Z*(H-1)/(H*(Z-1))*W,2*x+(2*H*Z-2*H-4*Z**2+2*Z)/Z*s)
images=(x,T.one(),s-x,T(W),(t*(s-x)-kap*es[2])/beta,x*W*t/beta+kap*(x*es[1]-es[0])/alpha,alpha/beta*W*t+kap*es[1],t,T.one());cache={};checks=[]
def ck(n,v):assert bool(v),n;checks.append(n)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def cv(c):
 if c not in cache:cache[c]=nk(ph(c.numerator())/ph(c.denominator()))
 return cache[c]
def nf(f):return T({tuple(map(int,ee)):c for ee,c in T(f).dict().items() if int(ee[0])+int(ee[1])<=1})
power={}
for j,im in enumerate(images):
 power[j,0]=T.one()
 for n in range(1,8):power[j,n]=nf(power[j,n-1]*im)
layers={}
for ee,c in d['G'].dict().items():
 if sum(ee)!=6:continue
 rd=sum(ee[6:]);kp=1-rd;f=T(cv(c))
 for j,n in enumerate(ee):
  if n:f=nf(f*power[j,int(n)])
 layers[kp]=layers.get(kp,T.zero())+f
answer=T.zero()
for kp,ff in sorted(layers.items()):
 if kp<0:
  ck('negative kappa cancellation '+str(kp),all(int(ee[3])>=-kp for ee in ff.dict()));ff=T({(int(ee[0]),int(ee[1]),int(ee[2]),int(ee[3])+kp):c for ee,c in ff.dict().items()})
 else:ff*=kap**kp
 answer+=ff
answer=T({tuple(map(int,ee)):nk(c) for ee,c in answer.dict().items() if c});print('FULL DUAL TERMS',len(answer.dict()),'sec',time.monotonic()-st,flush=True)
pure=T.hom([T.zero(),T.zero(),t,kap],T)(answer);sheet=T.hom([x,T.zero(),-2*kap,kap],T)(answer);first=T({(int(ee[0]),int(ee[1]),0,0):c for ee,c in answer.dict().items() if int(ee[2])==0 and int(ee[3])==1})
print('PURE',pure,'SHEET',sheet,'KAPPA FIRST T0',first,flush=True);ck('whole pure yw complete G source zero',pure==0);ck('whole lambda sheet complete G x dual zero',sheet==0);ck('whole G kappa first at t0 has no x/s linear',first==0)
U=b['ring'];Ub=U.base_ring();bh=Ub.ring().hom([k,r,H,Z],K);bg=U({tuple(map(int,ee)):Ub(c) for ee,c in b['Gbar'].dict().items()});gx=bg.derivative(U.gen(0));cf=sum((bh(c.numerator())/bh(c.denominator())*prod(v**int(n) for v,n in zip((T.zero(),T.zero(),T(W),T.zero()),ee)) for ee,c in gx.dict().items()),T.zero());ck('H cofactor constant at t0 vanishes',cf==0)
save({'ring':T,'dual_source':answer,'pure_remainder':pure,'lambda_sheet_remainder':sheet,'kappa_first_t0':first,'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_G30_J10_dual_reader.sobj'),compress=False);(p/'independent_subcritical_single_RR_G30_J10_dual_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS',len(checks),time.monotonic()-st,flush=True)
