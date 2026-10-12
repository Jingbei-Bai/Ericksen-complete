"""All-pitch exact whole-source Laurent leaders near B infinity and L zero.

Keeps a+b*U and all six pressure monomials separate until radical assembly.
Boundary coefficients alone are not asserted to be a no-lift certificate.
"""
from sage.all import *
from pathlib import Path
import time,resource,sys
OUT=Path(__file__).resolve().parent;START=time.time()
data=load(str(OUT/'trs_algebra_regular_frobenius.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
wall=int(sys.argv[sys.argv.index('--wall-order')+1]) if '--wall-order' in sys.argv else 0
weight=int(sys.argv[sys.argv.index('--L-weight')+1]) if '--L-weight' in sys.argv else 1
NF=QuadraticField(-1,'ii');ii=NF.gen()
R=PolynomialRing(NF,names=('w','p','j','P','I2','ell'));K=R.fraction_field()
w,p,j,P,I2,ell=map(K,R.gens());v=w*w
c=2-v*I2/j-P/j+v**3/j;d=v*P*I2/j**2-3*P/j;e=P**2/j**2
Q=PolynomialRing(K,'s');s=Q.gen();fpoly=1+c*s*s+d*s**4+e*s**6
PREC=16
def nk(z):
 z=K(z);n=z.numerator();den=z.denominator();g=n.gcd(den);return K(n//g)/K(den//g)
def cap(x,n):return Q([x[i] for i in range(n+1)])
if wall:
 FK=Q.fraction_field();ss=FK(s);bb=ss**-2
 aa=j/v-P/(v*bb);rr=I2-aa-j*(1+bb)/v;theta=-(j*j/v)*bb*FK(fpoly)
 ab=P/(v*bb**2);nn=v*(rr-aa)*ab-j*aa
 nf=4*theta*v*bb*ab**2+bb*nn**2+2*j*theta*(2*j/v-aa)
 dn=4*theta*ab+(rr-aa)*nn
 np=Q(-ss**2*nf/j);dp=Q(ss**2*dn)
 rat=[]
 for z in range(PREC):rat.append(nk((np[z]-sum(dp[i]*rat[z-i] for i in range(1,z+1)))/dp[0]))
 sqrtinv=sum(binomial(QQ(-1)/2,z)*(fpoly-1)**z for z in range(PREC//2+1))
 lrbase=cap(Q(rat)*sqrtinv,PREC-1)
 lrpoly=Q([nk(lrbase[z]) for z in range(PREC)])+ell/(-ii*w)*s**(wall-1)
 assert lrpoly[0]==1
 weight=1;Lunit=-ii*w
else:
 lrpoly=Q(1);Lunit=ell
cache={};pcache={}
def lpower(n,depth):
 key=(n,depth)
 if key not in pcache:
  z=Q(1)
  for _ in range(n):z=cap(z*cap(lrpoly,depth),depth)
  pcache[key]=z
 return pcache[key]
def small(a,t,l,depth):
 key=(a,t,l,depth)
 if key not in cache:
  ans=Q(0);power=Q(1)
  for z in range(depth//2+1):
   if z:power=cap(power*(fpoly-1),depth)
   ans+=binomial(QQ(t)/2,z)*power
  cache[key]=nk(cap(cap((1-P/j*s*s)**a,depth)*ans*lpower(l,depth),depth)[depth])
 return cache[key]
def leader(poly,caporder=16):
 if not poly:return None,K(0)
 dic=poly.dict();lowest=min(-m[4]-2*m[5]+weight*m[6] for m in dic)
 for order in range(lowest,lowest+caporder+1):
  ans=K(0)
  for mon,co in dic.items():
   ev,ep,ej,ea,et,eb,el=mon;base=-et-2*eb+weight*el
   if order<base:continue
   z=small(ea,et,el,order-base)
   if z:ans+=K(co)*v**ev*p**ep*j**ej*(j/v)**ea*(ii*j/w)**et*Lunit**el*z
  ans=nk(ans)
  if ans:return order,ans
 raise RuntimeError(('No leader',lowest,len(dic)))
results=[]
for row,C in enumerate(data['frobenius']):
 for odd,pol in enumerate(C):
  for pressure,f in sorted(pol.items()):
   if not f:continue
   on,cn=leader(f.numerator());od,cd=leader(f.denominator());coef=nk(cn/cd)
   item={'row':row,'U_parity':odd,'pressure':pressure,'order_without_U':on-od,'order_with_U':on-od-odd,'coefficient_without_U':coef}
   results.append(item)
   print(row,odd,pressure,'order',on-od-odd,'coefficient',coef.factor(),'seconds',round(time.time()-START,2),flush=True)
   save({'field':K,'wall_order':wall,'L_weight':weight,'results':results},str(OUT/('trs_B_infinity_L0_leaders_'+('wall'+str(wall) if wall else 'weight'+str(weight))+'.sobj')))
print('COMPLETE',round(time.time()-START,2),flush=True)
