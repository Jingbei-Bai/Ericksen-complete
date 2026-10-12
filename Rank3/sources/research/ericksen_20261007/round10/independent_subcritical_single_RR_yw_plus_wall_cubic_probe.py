"""True full-source cubic residual at the plus resonant yw wall."""
from sage.all import *
from pathlib import Path
import resource,os,time
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2)
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_full_regular_chart_H2.sobj'));g=load(str(p/'independent_subcritical_single_RR_full_complete_E_chart_H2.sobj'));P=d['ring'];OK=P.base_ring();K=PolynomialRing(QQ,names=('pc_r','pc_W')).fraction_field();r,W=K.gens();k=W**2;m=1+k;q=1+r;ph=OK.ring().hom([k,r],K);T=PolynomialRing(K,names=('pc_x','pc_z','pc_t','pc_kap','pc_u'));x,z,t,kap,u=T.gens();sources=d['sources']+[it['source'] for it in g['items']];cache={}
def cv(c):
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def nk(c):
 c=K(c)
 if not c:return K.zero()
 n,dd=c.numerator(),c.denominator();gg=n.gcd(dd);n,dd=n//gg,dd//gg;lc=dd.leading_coefficient();return K(n/lc)/K(dd/lc)
def nt(f,deg=3):return T({tuple(map(int,ee)):nk(c) for ee,c in T(f).dict().items() if sum(ee)<=deg and c})
def homogeneous(f,deg):return T({tuple(map(int,ee)):c for ee,c in f.dict().items() if sum(ee)==deg})
def mp(f,imgs):
 powers={};orders=[]
 for j,im in enumerate(imgs):
  orders.append(min((sum(ee) for ee in im.dict()),default=100));powers[j,0]=T.one()
  for n in range(1,16):powers[j,n]=nt(powers[j,n-1]*im)
 out=T.zero()
 for ee,c in f.dict().items():
  if sum(int(n)*o for n,o in zip(ee,orders))>3:continue
  val=T(cv(c))
  for j,n in enumerate(ee):
   if n:val=nt(val*powers[j,int(n)])
  out+=val
 return nt(out)
imgs=[x,z,W+u,(7*x+9*z)*(W+u),-(W+u)/2,16*x+14*z,27*r/q*(W+u)*(x+z),t,kap];initial=[mp(f,imgs) for f in sources];base=(K.zero(),K.zero(),W,K.zero(),-W/2,K.zero(),K.zero(),K.zero(),K.zero())
def em(f):return nk(sum((cv(c)*prod(v**int(n) for v,n in zip(base,ee)) for ee,c in f.dict().items() if not any(n and not v for n,v in zip(ee,base))),K.zero()))
A=matrix(K,len(sources),4,lambda i,j:em(sources[i].derivative(P.gen((3,4,5,6)[j]))));rows=(4,5,6,23);Ai=A.matrix_from_rows(rows).inverse();second=vector(T,[homogeneous(initial[i],2) for i in rows]);cor2=-Ai*second
for j,pos in enumerate((3,4,5,6)):imgs[pos]+=cor2[j]
corrected=[mp(f,imgs) for f in sources];range3=vector(T,[homogeneous(corrected[i],3) for i in rows]);residual=[nt(f-sum((T(A[i,j])*(second[j] if False else range3[j]) for j in range(4)),T.zero())) for i,f in enumerate(corrected)] if False else [nt(f-sum((T((A*Ai)[i,j])*range3[j] for j in range(4)),T.zero())) for i,f in enumerate(corrected)]
quadrics=[homogeneous(f,2) for f in residual];J=T.ideal([f for f in quadrics if f]);gb=J.groebner_basis();reduced=[homogeneous(f,3).reduce(gb) for f in residual]
print('QUAD GB',list(gb),'sec',time.monotonic()-st,flush=True)
for i,f in enumerate(reduced):
 if f:print('CUBIC',i,'terms',len(f.dict()),'FACT',str(f.factor())[:2500],flush=True)
save({'ring':T,'range_rows':rows,'range_matrix':A,'range2_correction':cor2,'images':imgs,'initial_jets':initial,'corrected_jets':corrected,'residuals':residual,'quadratic_ideal':J,'quadratic_GB':gb,'cubic_reductions':reduced,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_plus_wall_cubic_probe.sobj'),compress=False)
print('SAVED',time.monotonic()-st,flush=True)
