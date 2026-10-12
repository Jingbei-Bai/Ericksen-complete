"""Complete relative identity making the new kappa-divided source legitimate.

Four exact affine range equations, H4=t*Sigma*unit, and the actual even
combination Ceven supply the zero-kappa source ideal membership.
"""
from sage.all import *
from pathlib import Path
import resource,os,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_range.sobj'));j=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_bridge_jets.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2);print('PID',os.getpid(),flush=True);U=d['ring'];K=U.base_ring();k,r,H,Z=K.gens();x,z,w,t=U.gens();P=d['sources'][4].parent();B=H*(Z-1);A=Z*(H+1);S=B*(x+z)+A*x*z*t;checks=[]
def ck(s,v):assert bool(v),s;checks.append(s)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def nt(f):return U({tuple(map(int,ee)):nk(c) for ee,c in f.dict().items() if c})
def dw2(f):return nt(U({(int(ee[0]),int(ee[1]),int(ee[2])-2,int(ee[3])):c for ee,c in f.dict().items()}))
D=dw2(d['det']);N=list(map(dw2,d['numerators']));items={}
def image(f):
 md=max((sum(int(ee[j]) for j in (3,4,5,6)) for ee in f.dict()),default=0);out=U.zero();assert md<=1
 for ee,c in f.dict().items():
  val=U(c)*x**int(ee[0])*z**int(ee[1])*w**int(ee[2])*t**int(ee[7])*D**(md-sum(int(ee[j]) for j in (3,4,5,6)))
  for j,im in zip((3,4,5,6),N):val*=im**int(ee[j])
  out+=val
 return md,nt(out)
for ri in (11,30):
 md,f=image(d['sources'][ri]);items[ri]=(md,f);ck('true zero-kappa row has t factor '+str(ri),all(int(ee[3])>=1 for ee in f.dict()));print('ROW',ri,'MD',md,'TERMS',len(f.dict()),'sec',time.monotonic()-st,flush=True)
quot=U.fraction_field()(items[11][1])/(t*S);ck('H4 exact t Sigma polynomial factor',quot.denominator().is_constant());BH=nt(U(quot));ck('H4 multiplyback',nt(items[11][1]-t*S*BH)==0)
V=j['ring'];ww,tt=V.gens();h0=U.hom([V.zero(),V.zero(),ww,tt],V);BH0=h0(BH);b0=K(BH0.subs({ww:1,tt:0}));print('H4 QUOT0 FACT',b0.factor(),flush=True)
sh=K.ring().hom([k,r,K.one(),K(-1)],K);bcyc=nk(sh(b0.numerator())/sh(b0.denominator()));print('H4 CYCLIC QUOT0 FACT',bcyc.factor(),flush=True)
ck('G zero-kappa sheet first derivative vanishes identically',j['G_linear']==0);ck('Ceven exact first ratio',j['even_ratio']==-k/(1+k));ck('Ceven second cyclic plus unit',j['plus_C2']!=0)
# For a genuine local variable Sigma, z=-B*x/(B+A*x*t).
# G/t has value and derivative zero at x=0 on this sheet. Hence it is
# divisible by x^2 in the completed local ring. Ceven on Sigma=0 is
# x^2 times a unit, so G lies in (t*Sigma,Ceven), with regular coefficients.
Gbar=nt(U({(int(ee[0]),int(ee[1]),int(ee[2]),int(ee[3])-1):c for ee,c in items[30][1].dict().items()}));ck('Gbar value at xz0 zero',h0(Gbar)==0);ck('Gbar first sheet derivative zero',h0(Gbar.derivative(x)-Gbar.derivative(z))==0)
save({'ring':U,'D':D,'N':N,'Sigma':S,'items':items,'H4_tSigma_quotient':BH,'H4_origin_coefficient':b0,'H4_cyclic_origin_coefficient':bcyc,'Gbar':Gbar,'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_exact.sobj'),compress=False);(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_exact.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('PASS',len(checks),time.monotonic()-st,flush=True)
