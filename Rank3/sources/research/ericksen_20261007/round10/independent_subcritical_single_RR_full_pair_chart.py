"""Three complete genuine pair syzygies, regular on the full single-RR p13 chart."""
from sage.all import *
from pathlib import Path
import os,resource,time,sys,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));R=d['ring'];OK=R.base_ring();print('SOURCE FIELD',OK,flush=True);resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2);print('PID',os.getpid(),flush=True)
fixed=len(sys.argv)>1;K=PolynomialRing(QQ,names=('fp_k','fp_r') if fixed else ('fp_k','fp_r','fp_H','fp_Z')).fraction_field();k,r=K.gens()[:2];H,Z=(K(2),K(-1)) if fixed else K.gens()[2:];m=1+k;q=1+r
assert OK.ngens()==4
ph=OK.ring().hom([k/r,k*q/(r*m),H/r,Z*q/r],K);cache={}
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def cv(c):
 if c not in cache:cache[c]=nk(ph(c.numerator())/ph(c.denominator()))
 return cache[c]
P=PolynomialRing(K,names=('fp_x','fp_z','fp_w','fp_e0','fp_e1','fp_e2','fp_delta','fp_t','fp_kappa'));x,z,w,e0,e1,e2,de,t,kap=P.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));imgs=(x,P.one(),z,w,(t*z-kap*e2)/beta,x*w*t/beta+kap*(x*e1-e0)/alpha,alpha/beta*w*t+kap*e1,t,P.one());checks=[];items=[]
for it in d['items']:
 ans=P.zero();neg=P.zero()
 for ee,c in it['full'].dict().items():
  pd,rd=sum(ee[:6]),sum(ee[6:]);ke,dd=9-2*pd-3*rd,4-pd-rd;assert dd>=0 and ke>=-1,(pd,rd,ke,dd)
  val=cv(c)*de**dd*prod(v**int(n) for v,n in zip(imgs,ee))
  if ke==-1:neg+=val
  else:ans+=kap**ke*val
 assert all(int(ee[8])>=1 for ee in neg.dict()),it['pair'];ans+=P({tuple(int(nn)-(1 if j==8 else 0) for j,nn in enumerate(ee)):nk(c) for ee,c in neg.dict().items()});ans=P({tuple(map(int,ee)):nk(c) for ee,c in ans.dict().items()});checks.append('complete negative layer exactly divisible by kappa '+str(it['pair']));items.append({'pair':it['pair'],'source':ans,'negative_layer_before_division':neg,'original_source_combination':it['full_source_combination']});print('PAIR',it['pair'],'terms',len(ans.dict()),'sec',time.monotonic()-st,flush=True)
stem='independent_subcritical_single_RR_full_pair_chart'+('_H2' if fixed else '');save({'ring':P,'items':items,'chart_images':imgs,'normalizing_exponents':(9,4),'checks':checks,'scope':__doc__},str(p/(stem+'.sobj')),compress=False);(p/(stem+'.json')).write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n');print('PASS',time.monotonic()-st,flush=True)
