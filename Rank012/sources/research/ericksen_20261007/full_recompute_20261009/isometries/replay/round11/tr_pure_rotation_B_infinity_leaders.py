"""Exact Laurent leaders on finite-L B=infinity, full material constants.
This is a boundary calculation only, not an atlas completion.
"""
from sage.all import *
from pathlib import Path
import time,resource,sys
OUT=Path(__file__).resolve().parent;START=time.time()
LWEIGHT=int(sys.argv[sys.argv.index('--L-weight')+1]) if '--L-weight' in sys.argv else 0
NF=QuadraticField(-1,'ii');ii=NF.gen()
R=PolynomialRing(NF,names=('w','j','P','I2','L'));K=R.fraction_field()
w,j,P,I2,L=map(K,R.gens());v=w*w
c=2-v*I2/j-P/j+v**3/j;d=v*P*I2/j**2-3*P/j;e=P**2/j**2
Q=PolynomialRing(K,'s');s=Q.gen();ft=1+c*s+d*s*s+e*s**3
cache={}
def cap(p,n):return Q([p[i] for i in range(n+1)])
def small(a,t,delta):
 key=(a,t,delta)
 if key not in cache:
  total=Q(0);power=Q(1)
  for k in range(delta+1):
   if k:power=cap(power*(ft-1),delta)
   total+=binomial(QQ(t)/2,k)*power
  cache[key]=K(((1-P/j*s)**a*total)[delta])
 return cache[key]
def leader(poly,caporder=14):
 dic=poly.dict();lowest=min(-mon[3]-2*mon[4]+LWEIGHT*mon[5] for mon in dic)
 for order in range(lowest,lowest+caporder+1):
  ans=K(0)
  for mon,coeff in dic.items():
   base=-mon[3]-2*mon[4]+LWEIGHT*mon[5]
   if order<base or (order-base)%2:continue
   vv,jj,aa,tt,bb,ll=mon
   z=small(aa,tt,(order-base)//2)
   if z:ans+=K(coeff)*w**(2*vv)*j**jj*(j/v)**aa*(ii*j/w)**tt*L**ll*z
  if ans:
   return order,ans
 print('LEADER_NOT_FOUND',lowest,caporder,len(dic),flush=True)
 return None,None
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
names=['Bcompat_frobenius'] if '--Bcompat' in sys.argv else ['logU_frobenius','L_frobenius','flatness']
allout={}
for name in names:
 data=load(str(OUT/('tr_pure_rotation_'+name+'.sobj')))
 chosen=['1'] if '--constant-only' in sys.argv else list(data['coefficients'])
 subout={}
 for label in chosen:
  f=data['coefficients'][label];on,cn=leader(f.numerator());od,cd=leader(f.denominator())
  if cn is None or cd is None:continue
  extra=1 if label in ['EC','ED'] else 0
  coeff=cn/cd
  if '--details' in sys.argv: print('RAW',name,label,'numerator',on,cn.factor(),'denominator',od,cd.factor(),flush=True)
  subout[label]={'order_without_E':on-od,'order_with_E':on-od+extra,'coefficient_without_E':coeff}
  print(name,label,'order',on-od+extra,'coefficient_without_E',coeff.factor(), 'seconds',round(time.time()-START,2),flush=True)
 allout[name]=subout
 save({'field':K,'L_weight':LWEIGHT,'leaders':allout},str(OUT/('tr_pure_rotation_B_infinity_leaders'+('' if LWEIGHT==0 else '_Lweight'+str(LWEIGHT))+('_Bcompat' if '--Bcompat' in sys.argv else '')+'.sobj')))
print('COMPLETE',round(time.time()-START,2),flush=True)
