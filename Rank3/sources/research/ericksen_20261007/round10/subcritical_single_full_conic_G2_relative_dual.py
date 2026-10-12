"""Full Q6 source relative inclusion using the actual nonreduced conic.

Compute in C(H,Z,w)[x]/(x²-c)[s]/s², with the exact third-even
graph equation lifted by its derivative. Both coefficients must vanish.
"""
from sage.all import *
from pathlib import Path
import time,sys,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();g=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'));d=load(str(p/'independent_critical_single_RR3_normal.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));O=g['ring'];old=O.base_ring();oa,ob,oh,oz=old.gens();ov=O.gens();OP=d['ring'];OK=OP.base_ring();cyclic='--cyclic' in sys.argv
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
names=('ds_k','ds_r','ds_w') if cyclic else ('ds_k','ds_r','ds_H','ds_Z','ds_w');K=PolynomialRing(QQ,names=names).fraction_field()
if cyclic:k,r,w=K.gens();H=K(2);Z=K(-1)
else:k,r,H,Z,w=K.gens()
q=1+r;m=1+k;U=PolynomialRing(K,'ds_x');xx=U.gen();c=2*m/k*(k+w*w);E=U.quotient(xx*xx-c,'ds_b');x=E.gen();ph=old.ring().hom([k/r,k*q/(r*m),H/r,Z*q/r],K);po=OK.ring().hom([k/r,r,H,Z],K);cache={}
def nk(a):
 a=K(a)
 if not a:return K.zero()
 n,d0=a.numerator(),a.denominator();gg=n.gcd(d0);n,d0=n//gg,d0//gg;lc=d0.leading_coefficient();return K(n/lc)/K(d0/lc)
def ne(v):return E([nk(v[0]),nk(v[1])])
def ie(v):
 den=nk(v[0]*v[0]-c*v[1]*v[1]);assert den;return E([nk(v[0]/den),nk(-v[1]/den)])
def cv(a):
 a=old(a)
 if a not in cache:cache[a]=nk(ph(a.numerator())/ph(a.denominator()))
 return E(cache[a])
def oc(a):
 a=OK(a);return E(nk(po(a.numerator())/po(a.denominator())))
def mp0(f):
 ans=E.zero()
 for ex,a in f.dict().items():
  ex=tuple(map(int,ex));assert not any(ex[3:]);ans=ne(ans+oc(a)*x**ex[0]*(-x)**ex[1]*w**ex[2])
 return ans
state_ox,state_oz=OP.gens()[:2];C=d['graph_boundary'][27];Cz=mp0(C.derivative(state_oz));Cxminus=ne(mp0(C.derivative(state_ox))-Cz);xs=ne(-Cz*ie(Cxminus));print('LIFT xs',xs,'sec',time.monotonic()-st,flush=True)
def plus(a,b):return (ne(a[0]+b[0]),ne(a[1]+b[1]))
def times(a,b):return (ne(a[0]*b[0]),ne(a[0]*b[1]+a[1]*b[0]))
def sc(a,b):return (ne(a[0]*b),ne(a[1]*b))
one=(E.one(),E.zero());zero=(E.zero(),E.zero());pv=((x,xs),(-x,1-xs),(E(w),E.zero()))
def mpdual(f):
 ans=zero
 for ex,a in f.dict().items():
  ex=tuple(map(int,ex));assert not any(ex[3:]);v=one
  for image,n in zip(pv,ex[:3]):
   for unused in range(n):v=times(v,image)
  ans=plus(ans,sc(v,oc(a)))
 return ans
assert mpdual(C)==zero
ee=[mpdual(v) for v in d['graph_images'][3:6]];delta=mpdual(d['graph_images'][6]);alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));images=(pv[0],one,pv[1],pv[2],sc(ee[2],-1/beta),sc(plus(times(pv[0],ee[1]),sc(ee[0],-1)),1/alpha),ee[1],zero,one);kaporders=(0,0,0,0,1,1,1,0,0)
f=next(it['full_source'] for it in ev['items'] if int(it['Q_index'])==6);qv=f.parent().gens();it=g['items'][2];parts=[];keep=[]
for j in range(9):
 df=f.derivative(qv[j]);keep.append(df)
 if df:parts.append((it['field_numerators'][j],df))
for j,bp in enumerate((OK.gen(2),OK.gen(3))):
 scale=OK.gen(1) if j==0 else OK.gen(1)/(1+OK.gen(1));df=f.parent()({tuple(map(int,ex)):scale*a.derivative(bp) for ex,a in f.dict().items() if a.derivative(bp)});keep.append(df);parts.append((it['base_numerators'][j],df))
fc={}
def fmap(f):
 ky=id(f)
 if ky in fc:return fc[ky]
 out={}
 for ex,a in f.dict().items():
  ex=tuple(map(int,ex))
  if any(n and image==zero for n,image in zip(ex,images)):continue
  dg=sum(ex);rd=sum(ex[6:]);ki=sum(n*l for n,l in zip(ex,kaporders));v=one
  for image,n in zip(images,ex):
   for unused in range(n):v=times(v,image)
  ky=(dg,rd,ki);convert=cv if f.parent().base_ring()==old else oc;out[ky]=plus(out.get(ky,zero),sc(v,convert(a)))
 fc[id(f)]=out;return out
ans={}
for part,(left,right) in enumerate(parts):
 for (di,ji,ki),vi in fmap(left).items():
  for (dj,jj,kj),vj in fmap(right).items():
   dg=di+dj;rd=ji+jj;kp=13-2*dg-rd+ki+kj;dp=6-dg
   if kp>0:continue
   assert dp>=0;v=times(vi,vj)
   for unused in range(dp):v=times(v,delta)
   ans[kp]=plus(ans.get(kp,zero),v)
 print('PART',part,'sec',time.monotonic()-st,flush=True)
ans={kk:v for kk,v in ans.items() if v!=zero};print('REMAINDER',ans,'sec',time.monotonic()-st,flush=True)
stem='subcritical_single_full_conic_G2_relative_dual'+('_cyclic' if cyclic else '')
save({'ring':E,'conic_lift_xs':xs,'remainder':ans,'scope':__doc__,'cyclic':cyclic,'Q_source':6,'gauge_item':2},str(p/(stem+'.sobj')),compress=False)
print('DONE',time.monotonic()-st,flush=True)
