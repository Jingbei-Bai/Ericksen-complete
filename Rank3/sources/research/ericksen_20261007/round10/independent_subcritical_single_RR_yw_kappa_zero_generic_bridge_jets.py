"""Exact relative saturation bridge using x,z jets only, with w,t untruncated."""
from sage.all import *
from pathlib import Path
import resource,os,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_range.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+400*1024**2,)*2);print('PID',os.getpid(),flush=True);U=d['ring'];K=U.base_ring();k,r,H,Z=K.gens();m=1+k;q=1+r;x,z,w,t=U.gens();P=d['sources'][4].parent();px,pz,pw,pe0,pe1,pe2,pde,pt=P.gens();evs=(pe0,pe1,pe2,pde)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def dw2(f):
 assert all(int(ee[2])>=2 for ee in f.dict());return U({(int(ee[0]),int(ee[1]),int(ee[2])-2,int(ee[3])):c for ee,c in f.dict().items()})
D=dw2(d['det']);N=list(map(dw2,d['numerators']));B=H*(Z-1);A=Z*(H+1);S=B*(x+z)+A*x*z*t
V=PolynomialRing(K,names=('bj_w','bj_t'));ww,tt=V.gens();h0=U.hom([V.zero(),V.zero(),ww,tt],V)
def nv(f):return V({tuple(map(int,ee)):nk(c) for ee,c in V(f).dict().items() if c})
D0=nv(h0(D));base=(V.zero(),-2*Z*(H-1)/(H*(Z-1))*ww,V.zero(),V.zero());ef=P.hom([V.zero(),V.zero(),ww]+list(base)+[tt],V)
assert all(nv(h0(f)-a*D0)==0 for f,a in zip(N,base))
Dd=nv(h0(D.derivative(x)-D.derivative(z)));Ddd=nv(h0(D.derivative(x,2)+D.derivative(z,2)-2*D.derivative(x).derivative(z)));Dz=nv(h0(D.derivative(z)))
Ed=[nv(h0(f.derivative(x)-f.derivative(z))-a*Dd) for f,a in zip(N,base)];Ez=[nv(h0(f.derivative(z))-a*Dz) for f,a in zip(N,base)]
E2=[nv((h0(f.derivative(x,2)+f.derivative(z,2)-2*f.derivative(x).derivative(z))-a*Ddd)*D0/2-ee*Dd+A/B*tt*ez*D0) for f,a,ee,ez in zip(N,base,Ed,Ez)]
def residual_jets(f):
 assert ef(f)==0
 fd=f.derivative(px)-f.derivative(pz);fdd=fd.derivative(px)-fd.derivative(pz)
 lin=nv(ef(fd)*D0+sum((ef(f.derivative(v))*ee for v,ee in zip(evs,Ed)),V.zero()))
 quad=ef(fdd)*D0**2/2+A/B*tt*ef(f.derivative(pz))*D0**2
 quad+=sum((ef(fd.derivative(v))*ee*D0+ef(f.derivative(v))*e2 for v,ee,e2 in zip(evs,Ed,E2)),V.zero())
 quad+=sum((ef(f.derivative(v).derivative(vv))*ee*eex/2 for v,ee in zip(evs,Ed) for vv,eex in zip(evs,Ed)),V.zero())
 return lin,nv(quad)
result={}
for ri in (11,25,26,27,30):
 l,c=residual_jets(d['sources'][ri]);result[ri]=(l,c);print('ROW',ri,'LINEAR TERMS',len(l.dict()),'QUAD TERMS',len(c.dict()),'sec',time.monotonic()-st,flush=True)
ratio=V.fraction_field()(result[27][0])/result[25][0];print('RATIO',ratio.factor(),flush=True);assert ratio.denominator().is_constant();rat=K(ratio)
Ce2=nv(result[27][1]-rat*result[25][1]);LG=result[30][0];print('G LINEAR ZERO',LG==0,'G LINEAR',LG if LG else '',flush=True)
def plus(f):
 out=K.zero()
 for ee,c in V(f).dict().items():
  if ee[1]:continue
  assert int(ee[0])%2==0;out+=c*k**(int(ee[0])//2)
 return nk(out)
def spec(c,hh,zz):
 c=K(c);sh=K.ring().hom([k,r,K(hh),K(zz)],K);return nk(sh(c.numerator())/sh(c.denominator()))
c2plus=plus(Ce2);dplus=plus(D0)
for nm,val in [('D0',dplus),('Ce2',c2plus)]:
 print(nm,'FACT',val.factor(),flush=True);print('CYCLIC',nm,spec(val,1,-1).factor(),flush=True)
save({'ring':V,'D0':D0,'base':base,'graph_first':Ed,'graph_second':E2,'residual_jets':result,'even_ratio':rat,'Ceven_second_numerator':Ce2,'G_linear':LG,'plus_D0':dplus,'plus_C2':c2plus,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_kappa_zero_generic_bridge_jets.sobj'),compress=False);print('SAVED',time.monotonic()-st,flush=True)
