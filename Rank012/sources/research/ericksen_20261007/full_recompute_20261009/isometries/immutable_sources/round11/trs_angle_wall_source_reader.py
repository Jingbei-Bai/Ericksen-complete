"""Pull every actual archived first-order coefficient into the formal wall.
Independent from the hand-derived Schur/derivative formulas in the probe.
"""
from pathlib import Path
import os,resource,time,json
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from sage.all import *
root=Path(__file__).resolve().parent; start=time.time()
d=load(str(root/'trs_angle_wall_asymptotic_probe.sobj'))
source=load(str(root/'trs_algebra_first_order_shape.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+600*1024**2,vm+600*1024**2))
SR=d['M'].base_ring();K=SR.base_ring();e=SR.gen();w,p,j,I1,I2,KC,KB=map(K,K.ring().gens());v=w*w;q=v-p*p
PREC=d['PREC']; images=[SR(v),SR(p),SR(j),d['A'],d['T'],d['B'],d['L']]
checks=[]
def nk(x):
 x=K(x);n=x.numerator();de=x.denominator();g=n.gcd(de);n=n//g;de=de//g;lc=de.leading_coefficient()
 return K(n/lc)/K(de/lc)
def cl(x):
 x=SR(x);prec=min(PREC,x.precision_absolute())
 if not x:return SR(0).add_bigoh(prec)
 return sum(SR(nk(x[i]))*e**i for i in range(int(x.valuation()),int(prec)) if x[i]).add_bigoh(prec)
cache={}
def pw(i,n):
 key=(i,n)
 if key not in cache:cache[key]=cl(images[i]**n)
 return cache[key]
def poly(z):
 out=SR(0)
 for ex,c in z.dict().items():
  out += c*prod(pw(i,int(n)) for i,n in enumerate(ex) if n)
 return cl(out)
def pull(z):
 return cl(poly(z.numerator())/poly(z.denominator()))
def eq(name,x,y):
 z=cl(x-y)
 if z:raise AssertionError((name,z))
 checks.append(name)
 if len(checks)%10==0:print('checks',len(checks),'seconds',time.time()-start,flush=True)

eq('material I1',pull(source['I1']),SR(I1))
eq('material I2',pull(source['I2']),SR(I2))
eq('angle wall',pull(source['wall']),SR(0))
eq('U squared',pull(source['U_squared']),d['U2'])
M=matrix(SR,4,4)
for r in range(4):
 for c in range(4):
  M[r,c]=pull(source['M'][r,c]);eq(f'full M {r}{c}',M[r,c],d['M'][r,c])
rhs={}
for name,vec in source['rhs'].items():
 rhs[name]=vector(SR,[pull(x) for x in vec])
 for r in range(4):eq(f'full rhs {name} {r}',rhs[name][r],d['rhs_raw'][name][r])

# Independent 4x4 cofactor left null vector, column3, versus the Schur path.
# We use the first Schur compatibility because its leading coefficients are
# lower-order; all coefficients are recomputed from the pulled archive.
E=M[:2,:2];detE=cl(E.det())
eq('protected 2x2 determinant',detE,q*(v+d['L']**2)/v)
Einv=matrix(SR,[[E[1,1],-E[0,1]],[-E[1,0],E[0,0]]])/detE
N=(M[2:,2:]-M[2:,:2]*Einv*M[:2,2:]).apply_map(cl)
rr={key:(vec[2:]-M[2:,:2]*Einv*vec[:2]).apply_map(cl) for key,vec in rhs.items()}
c={key:cl(N[0,0]*vec[1]-N[1,0]*vec[0]) for key,vec in rr.items()}
ac=cl(c['constant']+KC*c['KC'])
dc=cl(KC*c['KC_over_U']+KB*c['KB_over_U']+c['one_over_U']+c['U']*d['U2'])
norm=cl(dc*dc-d['U2']*ac*ac)
assert norm.precision_absolute()>-2
assert nk(norm[-4]-j*j*q*KC*KC/v**5)==0
assert nk(norm[-3])==0
assert nk(norm[-2].subs({KC:0})+36*j*j*q*q/v**4)==0
checks.extend(['U-norm B2 coefficient','no B3/2 coefficient','KC0 U-norm B coefficient'])
out={'status':'PASS','checks':len(checks),'names':checks,'formal_variable':'B=e^-2, v=w^2, T~i*j/(w*e)',
 'norm_precision':str(norm.precision_absolute()),
 'leading_KC_nonzero':'[B^2] norm=j^2*(v-p^2)*KC^2/v^5',
 'leading_KC_zero':'[B] norm=-36*j^2*(v-p^2)^2/v^4',
 'scope':'source-algebra nonzero certificate; not a real large-B solution assumption'}
(root/'trs_angle_wall_source_reader.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
save({'source':str(root/'trs_algebra_first_order_shape.sobj'),'M':M,'rhs':rhs,'norm':norm,'c':c,'N':N,'checks':checks},str(root/'trs_angle_wall_source_reader.sobj'))
print('PASS',len(checks),'source and coefficient checks',flush=True)
