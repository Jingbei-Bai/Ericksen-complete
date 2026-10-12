"""All complete regular G normals on the cyclic four-root conic wall.
No no-lift conclusion until relative ideal inclusion is audited.
"""
from sage.all import *
from pathlib import Path
import os,resource,time
p=Path(__file__).resolve().parent;st=time.monotonic()
g=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
dd=load(str(p/'subcritical_single_cyclic_full_conic_wall_probe.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2)
E=dd['ring'];W=E.base_ring();K=W.base_ring();k,r=K.gens();w=E(W.gen());x=E.gen();q=1+r;m=1+k
old=g['ring'].base_ring();OK=ev['ring'].base_ring();oa,orr,H0,Z0=OK.gens();ph=old.hom([k/r,k*q/(r*m),1/r,-q/r],K);po=OK.hom([k/r,r,1,-1],K)
R=PolynomialRing(K,names=ev['ring'].variable_names());vs=R.gens();cache={}
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,den=c.numerator(),c.denominator();gg=nn.gcd(den);nn,den=nn//gg,den//gg;lc=den.leading_coefficient();return K(nn/lc)/K(den/lc)
def ne(f):return E([W([nk(f[i][j]) for j in range(2)]) for i in range(2)])
def cp(f):
 field=f.parent().base_ring();hh=ph if field==old else po;out={}
 for ee,c in f.dict().items():
  ky=(str(field),c)
  if ky not in cache:cache[ky]=hh(c)
  out[tuple(map(int,ee))]=cache[ky]
 return R(out)
def mul(a,b):
 out=[E.zero()]*4
 for i in range(4):
  for j in range(4-i):
   if a[i] and b[j]:out[i+j]=ne(out[i+j]+a[i]*b[j])
 return out
def add(a,b):return [ne(u+v) for u,v in zip(a,b)]
def sc(a,c):return [ne(z*c) for z in a]
one=[E.one()]+[E.zero()]*3;th=[E.zero(),E.one(),E.zero(),E.zero()]
ker=dd['items'][0]['kernel'][0]
directions=[vector(E,[0,0,0,0,0,0,-2,1]),vector(E,list(ker[:6])+[ker[6]-2,1])]
imageslist=[]
for v in directions:
 ee=[add(sc(one,z),sc(th,v[i])) for i,z in enumerate((-2*x*w,E.zero(),2*x))]
 tt=sc(th,v[6]);pv=[sc(one,x),one,sc(one,-x),sc(one,w),add(sc(tt,-x),sc(mul(th,ee[2]),-1)),add(sc(tt,x*w),mul(th,add(sc(ee[1],x),sc(ee[0],-1)))),add(sc(tt,w),mul(th,ee[1])),tt,one]
 imageslist.append(pv)
items=[]
for gi in (1,2,4):
 it=g['items'][gi];fv=list(map(cp,it['field_numerators']));bv=list(map(cp,it['base_numerators']))
 for qi in (1,5,6):
  OQ=next(t['full_source'] for t in ev['items'] if int(t['Q_index'])==qi);Q=cp(OQ)
  G=sum((fv[j]*Q.derivative(vs[j]) for j in range(9)),R.zero())
  for j,param in enumerate((H0,Z0)):
   scale=orr if j==0 else orr/(1+orr);df=OQ.parent()({tuple(map(int,e)):scale*c.derivative(param) for e,c in OQ.dict().items() if c.derivative(param)})
   G+=bv[j]*cp(df)
  values=[]
  for v,imgs in zip(directions,imageslist):
   ans={}
   for ee,c in G.dict().items():
    e=tuple(map(int,ee));dg=sum(e);rd=sum(e[6:]);dp=6-dg;assert dp>=0
    pow0=13-2*dg-rd+dp
    if pow0>1 or (dp and not v[3]):continue
    jet=one
    for im,n in zip(imgs,e):
     for unused in range(n):jet=mul(jet,im)
    for j,val in enumerate(jet):
     if pow0+j<=1 and val:ans[pow0+j]=ne(ans.get(pow0+j,E.zero())+E(c)*v[3]**dp*val)
   ans={j:z for j,z in ans.items() if z};assert all(j==1 for j in ans),(gi,qi,ans)
   values.append(ans.get(1,E.zero()))
  rec={'gauge':gi,'Q_index':qi,'normal_B':values[0],'normal_T':ne(values[1]-values[0]),'source':G}
  items.append(rec);print('G',gi,qi,'B',rec['normal_B'],'T',rec['normal_T'],'SECONDS',time.monotonic()-st,flush=True)
  save({'ring':E,'items':items,'directions':directions,'scope':__doc__},str(p/'subcritical_single_cyclic_conic_G_normal_cover.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
