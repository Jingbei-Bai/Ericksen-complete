"""True 3-jet of the exact kappa-saturated G30 source, after full four-range compensation."""
from sage.all import *
from pathlib import Path
import resource,os,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_yw_plus_wall_cubic_probe.sobj'));full=load(str(p/'independent_subcritical_single_RR_full_regular_chart_H2.sobj'));g=load(str(p/'independent_subcritical_single_RR_full_complete_E_chart_H2.sobj'));resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+450*1024**2,)*2);print('PID',os.getpid(),flush=True)
T=d['ring'];K=T.base_ring();r,W=K.gens();k=W**2;m=1+k;x,z,t,kap,u=T.gens();P=full['ring'];OK=P.base_ring();ph=OK.ring().hom([k,r],K);cache={};rows=d['range_rows'];A=d['range_matrix'];Ai=A.matrix_from_rows(rows).inverse();imgs=list(d['images']);sources=full['sources']+[it['source'] for it in g['items']]
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def cv(c):
 if c not in cache:cache[c]=nk(ph(c.numerator())/ph(c.denominator()))
 return cache[c]
def nt(f,n=4):return T({tuple(map(int,ee)):nk(c) for ee,c in f.dict().items() if sum(ee)<=n and c})
def hom(f,n):return T({tuple(map(int,ee)):c for ee,c in f.dict().items() if sum(ee)==n})
def mp(f):
 powers={};orders=[]
 for j,im in enumerate(imgs):
  orders.append(min((sum(ee) for ee in im.dict()),default=100));powers[j,0]=T.one()
  for nn in range(1,16):powers[j,nn]=nt(powers[j,nn-1]*im)
 out=T.zero()
 for ee,c in f.dict().items():
  if sum(int(n)*o for n,o in zip(ee,orders))>4:continue
  val=T(cv(c))
  for j,n in enumerate(ee):
   if n:val=nt(val*powers[j,int(n)])
  out+=val
 return nt(out)
cor3=-Ai*vector(T,[hom(d['corrected_jets'][i],3) for i in rows])
for j,pos in enumerate((3,4,5,6)):imgs[pos]+=cor3[j]
four=-Ai*vector(T,[hom(mp(sources[i]),4) for i in rows])
for j,pos in enumerate((3,4,5,6)):imgs[pos]+=four[j]
res=[mp(f) for f in sources];print('FOUR RANGE ZERO',all(res[i]==0 for i in rows),'sec',time.monotonic()-st,flush=True)
zero=T.hom([x,z,t,T.zero(),u],T);G0=zero(res[30]);H0=zero(res[11]);C=res[27]+k/m*res[25];C0=zero(C);a1=-5*W/m*(x-z);R=hom(G0-a1*H0,4)
sheet=T.hom([x,-x,t,T.zero(),u],T);Rs=sheet(R);assert all(int(ee[0])>=2 and int(ee[2])>=1 for ee in Rs.dict())
b1=T({(int(ee[0])-2,0,int(ee[2])-1,0,int(ee[4])):c*m/(16*k) for ee,c in Rs.dict().items()})
b2=t*b1;left=hom(R-b2*hom(C0,2),4);div=T.fraction_field()(left)/(t*(x+z)*(-2*k/5));assert div.denominator().is_constant();a2=T(div)
num=nt(res[30]-(a1+a2)*res[11]-b2*C);assert zero(num)==0;assert all(int(ee[3])>=1 for ee in num.dict())
Pstar=T({(int(ee[0]),int(ee[1]),int(ee[2]),int(ee[3])-1,int(ee[4])):c for ee,c in num.dict().items()})
print('A1',a1,'A2',a2,'B2',b2,'PSTAR',Pstar,'sec',time.monotonic()-st,flush=True)
save({'ring':T,'range_images4':imgs,'range4_correction':four,'residuals4':res,'cofactor_H_linear':a1,'cofactor_H_quadratic':a2,'cofactor_even_quadratic':b2,'saturated_source_jet3':Pstar,'scope':__doc__},str(p/'independent_subcritical_single_RR_yw_plus_wall_saturated_cubic_probe.sobj'),compress=False);print('SAVED',time.monotonic()-st,flush=True)
