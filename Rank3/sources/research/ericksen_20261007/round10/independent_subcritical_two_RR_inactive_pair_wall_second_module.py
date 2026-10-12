"""True complete paired wall module in seven WHOLE ideal generators.

Generators (mu*s,rho*s,kappa*s,rho*E,s²,sE,E²).  The module is evaluated
in a streaming quadratic material quotient, never a dense substituted rank.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();d=load(str(p/'independent_subcritical_two_RR_inactive_pair_second_module.sobj'));ki=load(str(p/'independent_subcritical_two_RR_inactive_pair_kappa_quadratic_module.sobj'));inc=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
F=d['ring'].base_ring();a,r,l=F.gens();q=1+r;m=1+a*r;C=PolynomialRing(F,names=('kap','rho'));k,rho=C.gens();CF=C.fraction_field();OK=inc['coefficient_field'];oph=OK.ring().hom([a,a*q/m,1/r,-q/r],F)
def ocv(c):return oph(c.numerator())/oph(c.denominator())
def iev(f):return sum((CF(ocv(co))*l**int(ex[0])*(rho/k)**int(ex[1])*(1/k)**int(ex[2]) for ex,co in f.dict().items()),CF.zero())
J=matrix(CF,6,6);J[0,0]=rho/l;J[0,3]=k/l;J[1,1]=1;J[2,0]=1;J[3,1]=rho;J[3,4]=-k;J[4,2]=1;J[5,5]=k*k
def crho(f):
 f=CF(f);nn=C(f.numerator());dd=C(f.denominator());dd0=dd.constant_coefficient();assert dd0,('regular original matrix',f)
 return F((nn.monomial_coefficient(rho)*dd0-nn.constant_coefficient()*dd.monomial_coefficient(rho))/dd0**2)
Nr=[]
for i,row in enumerate(inc['matrix'][:28].rows()):
 fac=k if i==0 or i>=10 else CF.one();nr=fac*vector(CF,[iev(f) for f in row])*J;Nr.append([crho(f) for f in nr])
pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pa['ring'];sx=S.gens();PK=S.base_ring();pph=PK.ring().hom([a,a*q/m,1/r,-q/r],F)
def pcv(c):return pph(c.numerator())/pph(c.denominator())
def ev(f):
 vals=[CF.one(),0,0,0,0,CF(l),l/k,rho/k,1/k];return sum((CF(pcv(co))*prod(CF(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),CF.zero())
for it in pa['items']:
 f=it['full'];top=S({ex:co for ex,co in f.dict().items() if sum(ex)==4});lower=S({ex:co for ex,co in f.dict().items() if sum(ex)==3});mr=vector(CF,[ev(top.derivative(sx[j])) for j in [1,2,3,4,6]]+[ev(lower)]);Nr.append([crho(z) for z in k*mr*J])
Nr=matrix(F,Nr);N=d['first_matrix'];Nk=ki['kappa_first_matrix'];nu=vector(F,[1,-1/l,0,0,0,0]);bfree=load(str(p/'independent_subcritical_two_RR_inactive_pair_free_boundary.sobj'));BP=bfree['ring'];ss,vv,EE=BP.gens();pure=[]
for f in bfree['rows']:
 f=BP(f.subs({vv:-ss/l}));pure.append([f.monomial_coefficient(mon) for mon in [ss**2,ss*EE,EE**2]])
lin=N*nu;mu=vector(F,[z.derivative(l)/(2*m*l) for z in lin]);rs=Nr*nu;ks=Nk*nu;re=Nr.column(2);raw=matrix(F,[[mu[j],rs[j],ks[j],re[j]]+pure[j] for j in range(31)])
print('TRUE original coefficient blocks ready',time.time()-st,flush=True)
B=PolynomialRing(QQ,names=('a','r')).fraction_field();aa,rr=B.gens();U=PolynomialRing(B,'l');ll=U.gen();ph=F.ring().hom([aa,rr,ll],U);items=[];checks=[]
for eps,erow in [(1,1),(-1,2)]:
 rranges=[0,erow,7,28];cols=[1,3,4,5];BN=N.matrix_from_columns(cols);A=BN.matrix_from_rows(rranges);comp=-A.inverse()*raw.matrix_from_rows(rranges);M=raw+BN*comp;c=eps*aa*rr/(1+aa*rr);wall=ll**2-c
 def cv(f):
  nn=ph(f.numerator())%wall;dd=ph(f.denominator())%wall;zz=(nn*dd.inverse_mod(wall))%wall;return (B(zz[0]),B(zz[1]))
 def add(x,y):return (x[0]+y[0],x[1]+y[1])
 def neg(x):return (-x[0],-x[1])
 def mul(x,y):return (x[0]*y[0]+c*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
 def inv(x):
  no=x[0]*x[0]-c*x[1]*x[1];assert no;return (x[0]/no,-x[1]/no)
 MR=[[cv(z) for z in row] for row in M.rows()];pivs=[];pcs=[];prs=[]
 for j,row in enumerate(MR):
  v=list(row)
  for pc,pv in zip(pcs,pivs):
   z=v[pc]
   if z!=(0,0):v=[add(x,neg(mul(z,y))) for x,y in zip(v,pv)]
  co=next((k for k,z in enumerate(v) if z!=(0,0)),None)
  if co is not None:
   iv=inv(v[co]);v=[mul(iv,z) for z in v];pcs.append(co);prs.append(j);pivs.append(v)
 print('WALL MODULE',eps,'rank',len(pcs),'cols',pcs,'rows',prs,flush=True)
 for j in [2 if eps==1 else 1,3,5,11,13,17,19,23,25,29]:
  print('WALL ROW',eps,j,[(str(x.factor()) if x else '0',str(y.factor()) if y else '0') for x,y in MR[j]],flush=True)
 items.append({'epsilon':eps,'range_rows':rranges,'range_matrix':A,'true_range_correction':comp,'full_module_before_wall':M,'material_wall':wall,'module_pairs':MR,'pivot_rows':prs,'pivot_columns':pcs,'reduced_pivot_rows':pivs})
save({'field':F,'original_first':N,'original_rho_first':Nr,'original_kappa_first':Nk,'whole_generators':['mu*s','rho*s','kappa*s','rho*E','s²','sE','E²'],'raw_source_module':raw,'items':items,'scope':'True complete cyclic all-material second module on both paired walls; whole ideal proof pending'},str(p/'independent_subcritical_two_RR_inactive_pair_wall_second_module.sobj'),compress=False)
print('END',time.time()-st,flush=True)
