"""Exact paired missing-RR second coefficients, with all l retained.
Only pure-kernel/tangent coefficients are evaluated; no normal Hessian surrogate.
"""
from sage.all import *
from pathlib import Path
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.time();inc=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
print('LinuxPID',os.getpid(),flush=True)
F=PolynomialRing(QQ,names=('a','r','l')).fraction_field();a,r,l=F.gens();q=1+r;m=1+a*r
ok=inc['coefficient_field'];oph=ok.ring().hom([a,a*q/m,1/r,-q/r],F)
def ocv(c):return oph(c.numerator())/oph(c.denominator())
C=PolynomialRing(F,names=('rho','kap'));rho,kap=C.gens();CF=C.fraction_field()
J=matrix(CF,6,6);J[0,0]=rho/l;J[0,3]=kap/l;J[1,1]=1;J[2,0]=1;J[3,1]=rho;J[3,4]=-kap;J[4,2]=1;J[5,5]=kap*kap
def iev(f):return sum((CF(ocv(co))*l**int(ex[0])*rho**int(ex[1])*kap**(-int(ex[1])-int(ex[2])) for ex,co in f.dict().items()),CF.zero())
def lim(f):
 f=CF(f);nn=f.numerator();dd=f.denominator();dd0=dd.constant_coefficient();assert dd0,('pole',f)
 return F(nn.constant_coefficient()/dd0)
N=[]
for i,row in enumerate(inc['matrix'][:28].rows()):
 fac=kap if i==0 or i>=10 else CF.one();rr=fac*vector(CF,[iev(f) for f in row])*J;N.append([lim(f) for f in rr])
P=PolynomialRing(F,names=('rho','kap','E0'));rho,kap,E0=P.gens();checks=[]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def tr(f,d=2):return P({ee:co for ee,co in P(f).dict().items() if sum(ee)<=d})
pv=[P.one(),P.zero(),P.zero(),P.zero(),P.zero(),P(l)];Rv=[l+kap*E0,rho,P.one()];vals=pv+Rv;es=[E0,P.zero(),P.zero()]
qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));QO=qp['ring'];QK=QO.base_ring();ph=QK.ring().hom([a,r,1,-1],F)
def qcv(c):return ph(c.numerator())/ph(c.denominator())
def pull(dd,cv,d=2):
 out=P.zero()
 for ee,co in dd.items():
  e=tuple(map(int,ee));e=e+(0,)*(9-len(e));mon=P.one()
  for j,n in enumerate(e):
   if n:mon=tr(mon*tr(vals[j]**n,d),d)
   if not mon:break
  if mon:out+=P(cv(co))*mon
 return tr(out,d)
def face(f,ij):return {tuple(map(int,ee)):co for ee,co in f.dict().items() if (sum(ee[:6]),sum(ee[6:]))==ij}
off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset'];DC=matrix(QO,qp['D'])*matrix(QO,qp['C']);Dd=matrix(QO,qp['D'])*vector(QO,off);rows=[]
Q=qp['full_Q0'];rows.append(tr(pull(qp['Q0_L_cofactors'][0].dict(),qcv)*E0+pull(face(Q,(2,1)),qcv)+kap*pull(face(Q,(3,0)),qcv)))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for it in ev['items']:
 rows.append(tr(P(qcv(it['square_unit']))*es[int(it['square_component'])]**2+pull(it['A_p_quadratic'][0].dict(),qcv)*E0+pull(it['pure_p4'].dict(),qcv)))
for j,f in enumerate(qp['full_R_sources']):rows.append(tr(pull(DC[j,0].dict(),qcv)*E0+pull(Dd[j].dict(),qcv)))
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
for j,it in enumerate(hi['items']):
 ff=tr(pull(it['L_cofactors_original'][0],ocv)*E0+pull(it['H41_original'],ocv))
 rows.append(ff);print('H',j,'terms',len(ff.dict()),'sec',round(time.time()-st,2),flush=True)
pair=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pair['ring'];Sx=S.gens();pk=S.base_ring();pkph=pk.ring().hom([a,a*q/m,1/r,-q/r],F)
def pcv(c):return pkph(c.numerator())/pkph(c.denominator())
def pf(f,ij):return S(face(f,ij))
def ev0(f):
 vv=[F.one(),0,0,0,0,l,l,0,1];return sum((pcv(co)*prod(F(vv[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),F.zero())
for j,it in enumerate(pair['items']):
 f=it['full'];f22=pf(f,(2,2));f31=pf(f,(3,1));f03=pf(f,(0,3))
 # New regular chart first coefficients (s=p23,v=p21,E0,e1,e2,delta).
 N.append([ev0(f31.derivative(Sx[3])),ev0(f31.derivative(Sx[2])),ev0(f22.derivative(Sx[6])),ev0(f22.derivative(Sx[1]))/l,-ev0(f22.derivative(Sx[4])),ev0(f03)])
 lead=pull(f22.dict(),pcv,3);ck('pair22 divisible kap '+str(j),all(int(ee[1])>=1 for ee in lead.dict()));rows.append(tr(lead/kap+pull(f31.dict(),pcv)+kap*pull(pf(f,(4,0)).dict(),pcv)))
N=matrix(F,N);ck('all31 first E0 zero',not N.column(2));range_cols=[0,1,3,4,5];NR=N.matrix_from_columns(range_cols)
for j,f in enumerate(rows):
 ck('complete pure-kernel two-generator span '+str(j),f==f.monomial_coefficient(rho*E0)*rho*E0+f.monomial_coefficient(E0*E0)*E0*E0)
items=[]
rr=(0,1,2,5,7);A=NR.matrix_from_rows(rr);print('RANGE',rr,A.det().factor(),flush=True)
comp=-A.inverse()*vector(P,[rows[j] for j in rr]);res=[rows[j]+sum((NR[j,h]*comp[h] for h in range(5)),P.zero()) for j in range(len(rows))];M=matrix(F,[[f.monomial_coefficient(rho*E0),f.monomial_coefficient(E0*E0)] for f in res]);print('RANK',M.rank(),flush=True)
rp=M.transpose().pivots();print('MODULE_PIVOTS',rp,flush=True)
if len(rp)==2:print('TWO_MINOR',M.matrix_from_rows(rp).det().factor(),flush=True)
print('DELTA_COMP',comp[4],flush=True)
for j,f in enumerate(res):
 if f:print('RES',j,f,flush=True)
save({'ring':P,'rows':rows,'first_matrix':N,'range_rows':rr,'range_matrix':A,'range_compensation':comp,'residuals':res,'module_matrix':M,'checks':checks,'scope':'Exact cyclic all-material paired second coefficients; relative source membership/generic-base lift pending'},str(p/'independent_subcritical_two_RR_inactive_pair_second_module.sobj'),compress=False)
print('PASS',len(checks),'sec',time.time()-st,flush=True)
