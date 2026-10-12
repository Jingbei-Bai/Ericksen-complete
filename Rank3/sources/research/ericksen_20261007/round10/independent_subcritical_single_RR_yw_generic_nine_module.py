"""Full generic H/Z nine-module and sparse determinant on the yw family."""
from sage.all import *
from pathlib import Path
import resource,os,time
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+550*1024**2,)*2);p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_single_RR3_normal.sobj'));OP=d['ring'];OK=OP.base_ring();K=PolynomialRing(QQ,names=('gn_k','gn_r','gn_H','gn_Z','gn_W')).fraction_field();k,r,H,Z,W=K.gens();ph=OK.ring().hom([k/r,r,H,Z],K);T=PolynomialRing(K,names=('gn_x','gn_z'));x,z=T.gens();cache={}
def nk(c):
 c=K(c)
 if not c:return K.zero()
 n,dd=c.numerator(),c.denominator();gg=n.gcd(dd);n,dd=n//gg,dd//gg;lc=dd.leading_coefficient();return K(n/lc)/K(dd/lc)
def cv(c):
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def mp(f,imgs):return sum((cv(c)*prod(v**int(n) for v,n in zip(imgs,ee)) for ee,c in f.dict().items() if not any(n and not v for n,v in zip(ee,imgs))),T.zero())
graph=[mp(f,(x,z,T(W),T.zero(),T.zero(),T.zero(),T.zero())) for f in d['graph_images']];base=tuple(f(0,0) for f in graph);A=matrix(K,28,4,lambda i,j:nk(mp(d['normal_matrix'][i,(3,4,5,6)[j]],base).constant_coefficient()));Q=matrix(K,28,7)
for i in range(28):
 for pc,col in ((8,0),(7,2)):
  f=mp(d['normal_matrix'][i,pc],graph);Q[i,col]=nk(f.monomial_coefficient(x));Q[i,col+1]=nk(f.monomial_coefficient(z))
 f=mp(d['graph_boundary'][i],(x,z,T(W),T.zero(),T.zero(),T.zero(),T.zero()))
 for j,mv in enumerate((x*x,x*z,z*z),4):Q[i,j]=nk(f.monomial_coefficient(mv))
 print('ROW',i,'sec',time.monotonic()-st,flush=True)
c=(3*k-W*W)/(k+W*W);B=matrix(K,5,7,[[1,1,0,0,0,0,0],[0,0,1,1,0,0,0],[0,0,0,0,1,2,1],[2,0,1,0,0,0,0],[0,0,0,0,1-c,-c,0]]);sel=matrix(K,7,5);sel[0,0]=1;sel[1,1]=1;sel[2,2]=1;sel[4,3]=1;sel[5,3]=-1;sel[6,4]=1;N=Q*sel*(B*sel).inverse();assert N*B==Q;M=A.augment(N);arows=(4,5,6,23);A0=A.matrix_from_rows(arows);det4=nk(A0.det());print('DET4',det4.factor(),'sec',time.monotonic()-st,flush=True);S=N-A*A0.inverse()*N.matrix_from_rows(arows);srows=(11,12,25,26,27);SS=S.matrix_from_rows(srows);assert all(not SS[i,j] for i in range(2) for j in range(2,5));dodd=nk(SS[0,0]*SS[1,1]-SS[0,1]*SS[1,0]);print('DET_ODD',dodd.factor(),'sec',time.monotonic()-st,flush=True);dev=nk(SS.matrix_from_rows_and_columns((2,3,4),(2,3,4)).det());print('DET_EVEN',dev.factor(),'sec',time.monotonic()-st,flush=True);det9=nk(det4*dodd*dev);print('DET9',det9.factor(),'sec',time.monotonic()-st,flush=True);save({'ring':T,'graph':graph,'range_matrix':A,'quadratic_matrix':Q,'basis':B,'selector':sel,'matrix':M,'range_rows':arows,'residual_rows':srows,'range_det':det4,'odd_det':dodd,'even_det':dev,'minor':det9,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_generic_nine_module.sobj'),compress=False)
