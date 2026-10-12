"""Raw hatFi(Q0) third jets WITH the true second range correction.

The calculation is over the two quadratic relations l^2=dd and P(y)=0.
All first normal/tangent compatibility coefficients are also retained, so
zero tangent corrections are never silently used as a no-lift argument.
"""
from sage.all import *
from pathlib import Path
import time,os,resource,json
p=Path(__file__).resolve().parent;st=time.monotonic();c=load(str(p/'independent_ratio_incoming12_nullwall_second_compatibility.sobj'));raw=load(str(p/'independent_ratio_raw_Q0_normal_two.sobj'));U=c['ring'];K=U.base_ring();a,b,h,z=K.gens();l=U.gen();dd=-c['relation'][0];P=c['primitive_quadratic'];rel0=-P[0]/P[2];rel1=-P[1]/P[2];r=(a-b)/(a*(b-1));q=1+r;Z=r*z/q;kap=(Z-1)/Z
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+700*1024**2,)*2)
# Basis (1,l,y,l*y), using raw K coefficients rather than nested fields.
zero=(K(0),)*4;one=(K(1),K(0),K(0),K(0));ell=(K(0),K(1),K(0),K(0));yy=(K(0),K(0),K(1),K(0))
def add(x,y):return tuple(x[i]+y[i] for i in range(4))
def scale(c,x):return tuple(K(c)*f for f in x)
def mul(x,y):
 out=[K(0)]*4
 for i,ci in enumerate(x):
  if not ci:continue
  for j,cj in enumerate(y):
   if not cj:continue
   lp=i%2+j%2;yp=i//2+j//2;coef=ci*cj*(dd if lp>=2 else 1);lp%=2
   if yp<2:out[lp+2*yp]+=coef
   else:out[lp]+=coef*rel0;out[lp+2]+=coef*rel1
 return tuple(out)
def powe(x,n):
 out=one
 while n:
  if n%2:out=mul(out,x)
  n//=2
  if n:x=mul(x,x)
 return out
def uv(f):
 f=U(f)%c['relation'];return (K(f[0]),K(f[1]),K(0),K(0))
def fromY(f):
 out=zero
 for i,cc in enumerate(f):out=add(out,mul(uv(cc),powe(yy,i)))
 return out
delta=U(c['range_denominator']);nn=delta[0]**2-dd*delta[1]**2;di=(delta[0]/nn,-delta[1]/nn,K(0),K(0));assert mul(uv(delta),di)==one
n2=[mul(fromY(f),di) for f in c['range_second_numerators']]
assert n2[3:]==[zero,zero,zero]
print('range sizes',[[len(f.numerator().dict()) if f else 0 for f in v] for v in n2],flush=True)
def smul(x,y):
 out={}
 for i,a0 in x.items():
  for j,b0 in y.items():
   if i+j<=3:out[i+j]=add(out.get(i+j,zero),mul(a0,b0))
 return {i:v for i,v in out.items() if any(v)}
def spow(x,n):
 out={0:one}
 for _ in range(n):out=smul(out,x)
 return out
k1=mul(ell,yy);p13=add(scale(kap,ell),scale(-1,k1))
states=[{0:one},{1:p13,2:n2[0]},{2:n2[1]},{2:n2[2]},{1:k1},{0:ell},{0:scale(-2,ell)},{},{0:scale(-2,one)}]
cache={};moncache={}
def factor(ref):
 if ref in cache:return cache[ref]
 out={}
 for ee,cc in raw['polynomial_pool'][ref].items():
  e=tuple(map(int,ee))
  if e[7]:continue
  if e not in moncache:
   val={0:one}
   for j,ej in enumerate(e):
    if ej:val=smul(val,spow(states[j],ej))
    if not val:break
   moncache[e]=val
  for ep,v in moncache[e].items():
   key=(sum(e),ep);out[key]=add(out.get(key,zero),scale(cc,v))
 cache[ref]={key:v for key,v in out.items() if any(v)};return cache[ref]
def product(f,g):
 out=zero
 for (d,i),v in f.items():
  for (e,j),w in g.items():
   tp=5-d-e
   if tp>=0 and i+j+tp==3:out=add(out,mul(v,w))
 return out
F=raw['ring'];lf,uf,vf=F.gens()
def fieldval(f):
 out=zero
 for ee,cc in F(f).dict().items():
  e=tuple(map(int,ee))
  if e[1]:continue
  out=add(out,scale(K(cc)*(-2)**e[2],powe(ell,e[0])))
 return out
n1=[p13,zero,zero,k1,zero,one];items=[]
for it in raw['items']:
 idx=int(it['direction']);co=it['coefficients'];first=[fieldval(f) for f in co[:6]];assert all(v==zero for v in first)
 tangent={}
 for name,tg in [('dl',lf),('du',uf),('dv',vf)]:
  val=zero
  for j in range(6):val=add(val,mul(fieldval(co[j].derivative(tg)),n1[j]))
  tangent[name]=val
 third=zero
 for left,right in it['full_product_DAG']:third=add(third,product(factor(left),factor(right)))
 items.append({'direction':idx,'third_with_range':third,'tangent_second_coefficients':tangent})
 print('source',idx,'third-nonzero',[j for j,v in enumerate(third) if v],'tangents',{key:any(v) for key,v in tangent.items()},'seconds',round(time.monotonic()-st,2),flush=True)
 save({'coefficient_field':K,'l_square':dd,'y_quadratic':P,'basis':['1','l','y','l*y'],'second_range':n2,'items':items,'scope':__doc__},str(p/'independent_ratio_incoming12_raw_third_range.sobj'),compress=False)
(p/'independent_ratio_incoming12_raw_third_range.json').write_text(json.dumps({'source_count':len(items),'nonzero_thirds':[[j for j,v in enumerate(it['third_with_range']) if v] for it in items],'tangent_nonzero':[{key:bool(any(v)) for key,v in it['tangent_second_coefficients'].items()} for it in items],'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n');print('DONE',time.monotonic()-st,flush=True)
