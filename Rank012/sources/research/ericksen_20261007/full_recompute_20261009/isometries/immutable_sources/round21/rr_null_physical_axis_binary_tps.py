"""Complete the square in chi and compute a short exact even-power expansion.

Every coefficient operation is canonically normalized. The output is not a
whole theorem until all leading factors and the transformed-coordinate guards
are reviewed.
"""
from sage.all import *
from pathlib import Path
import time,signal,json,sys,resource
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(140)
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1400*1024**2,vm+1400*1024**2))
tag=sys.argv[1] if len(sys.argv)>1 else 'plus';order=int(sys.argv[2]) if len(sys.argv)>2 else 3
d=load(str(O/f'rr_null_physical_axis_reduction_checkpoint_{tag}.sobj'));F=d['field'];P=F.ring();E=d['e_ring'];e=E.gen();N=d['reduced_matrix']
R=PolynomialRing(F,names=('a','b'));a,b=R.gens()
def cf(z):
 z=F(z)
 if not z:return z
 lc=z.denominator().leading_coefficient();return F(z.numerator()/lc)/F(z.denominator()/lc)
def norm(z):return R({m:cf(c) for m,c in R(z).dict().items()})
def val(z):return z.numerator().valuation()-z.denominator().valuation() if z else Infinity
def ec(z,n):
 if not z or val(z)>n:return F(0)
 nu=z.numerator();de=z.denominator();vn=nu.valuation();vd=de.valuation();vv=vn-vd;out=[]
 for j in range(n-vv+1):out.append(cf((nu[vn+j]-sum(de[vd+k]*out[j-k] for k in range(1,j+1)))/de[vd]))
 return out[-1]
vals=[min(val(z) for z in row) for row in N.rows()];N=diagonal_matrix(N.base_ring(),[e**(-v) for v in vals])*N
qq=[[[R(ec(z,2*j)) for j in range(order)] for z in row] for row in N.rows()]
assert all(ec(z,1)==0 and ec(z,3)==0 for z in N.list())
print('coefficient truncation',time.time()-st,flush=True)
zero=[R(0)]*order
def add(u,v):return [norm(x+y) for x,y in zip(u,v)]
def neg(u):return [-x for x in u]
def sub(u,v):return add(u,neg(v))
def scale(u,z):return [norm(x*z) for x in u]
def mul(u,v):
 out=[]
 for n in range(order):
  out.append(R(0))
  for j in range(n+1):out[n]=norm(out[n]+norm(u[j]*v[n-j]))
 return out
def inv(u):
 out=[norm(1/u[0])]
 for n in range(1,order):
  vv=R(0)
  for j in range(1,n+1):vv=norm(vv+norm(u[j]*out[n-j]))
  out.append(norm(-vv/u[0]))
 return out
G=qq[1][5];iG=inv(G);u=scale(mul(qq[1][2],iG),QQ(1)/2);v=scale(mul(qq[1][4],iG),QQ(1)/2)
uv=mul(u,v);uu=mul(u,u);vv=mul(v,v)
def newrow(row):
 A,B,C,D,Ee,H=row
 aa=sub(add(A,mul(H,uu)),mul(C,u))
 ab=sub(sub(add(B,scale(mul(H,uv),2)),mul(C,v)),mul(Ee,u))
 bb=sub(add(D,mul(H,vv)),mul(Ee,v))
 con=add(add(scale(aa,a*a),scale(ab,a*b)),scale(bb,b*b))
 lin=add(scale(sub(C,scale(mul(H,u),2)),a),scale(sub(Ee,scale(mul(H,v),2)),b))
 return con,lin
AA,MM=newrow(qq[0]);BB,ll=newrow(qq[1]);CC,NN=newrow(qq[2])
assert not any(ll)
CC=sub(CC,mul(mul(qq[2][5],iG),BB))
DD=sub(mul(MM,CC),mul(NN,AA));HH=add(mul(G,mul(AA,AA)),mul(BB,mul(MM,MM)))
def first(ss):
 for j,z in enumerate(ss):
  if z:return j,z
 return None,R(0)
orders={k:first(z)[0] for k,z in [('A',AA),('B',BB),('C',CC),('D',DD),('H',HH)]}
di,df=first(DD);hi,hf=first(HH)
save({'field':F,'ring':R,'order':order,'e_squared_orders':orders,'A':AA,'B':BB,'C':CC,'M':MM,'N':NN,'D':DD,'H':HH,'G':G,'u':u,'v':v},str(O/f'rr_null_physical_axis_binary_tps_{tag}.sobj'))
print('orders',orders,time.time()-st,flush=True)
if di is None or hi is None:
 print('MORE_PRECISION_REQUIRED',flush=True);sys.exit(0)
def primitive(f):
 den=lcm(z.denominator() for z in f.coefficients());nn=[P(z*den) for z in f.coefficients()];gg=gcd(nn);sc=cf(F(gg)/F(den));return R(f/sc),sc
dp,ds=primitive(df);hp,hs=primitive(hf)
flat=PolynomialRing(P.base_ring(),names=P.variable_names()+('a','b'));ss,tt,av,bv=flat.gens();bf=P.hom([ss,tt],flat)
# Use explicit polynomial coefficient conversion rather than fraction-field hom.
def flatten(f):
 return sum(bf(P(co))*av**mon[0]*bv**mon[1] for mon,co in f.dict().items())
fd=flatten(dp).factor();fh=flatten(hp).factor()
print('D factors',fd,flush=True);print('H factors',fh,flush=True)
U=PolynomialRing(P,'a');fv=U([P(hp.monomial_coefficient(a**j*b**(4-j))) for j in range(5)]);gv=U([P(dp.monomial_coefficient(a**j*b**(3-j))) for j in range(4)])
res=fv.resultant(gv);fr=res.factor() if res else 0
save({'resultant':res,'factorization':fr,'D':dp,'H':hp,'D_factorization':fd,'H_factorization':fh,'D_scale':ds,'H_scale':hs},str(O/f'rr_null_physical_axis_binary_resultant_{tag}.sobj'))
rep={'status':'EXACT_BINARY_INITIAL_RESULTANT_REQUIRES_GUARD_AUDIT','tag':tag,'orders':orders,'resultant':str(fr),'D_factorization':str(fd),'H_factorization':str(fh),'D_scale':str(ds.factor()),'H_scale':str(hs.factor()),'seconds':time.time()-st}
(O/f'rr_null_physical_axis_binary_resultant_{tag}.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
