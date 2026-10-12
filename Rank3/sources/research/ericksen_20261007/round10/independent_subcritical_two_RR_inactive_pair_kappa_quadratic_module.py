"""True kappa*E0^2 coefficients using complete original incoming two-jets.
Includes kappa variation of every range row and true range compensation.
"""
from sage.all import *
from pathlib import Path
import os,resource,time
p=Path(__file__).resolve().parent;st=time.time();d=load(str(p/'independent_subcritical_two_RR_inactive_pair_second_module.sobj'));inc=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));two=load(str(p/'independent_ratio_incoming12_generic_normal_two.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
F=d['ring'].base_ring();a,r,l=F.gens();q=1+r;m=1+a*r;C=PolynomialRing(F,'kap');k=C.gen();CF=C.fraction_field();OK=inc['coefficient_field'];oph=OK.ring().hom([a,a*q/m,1/r,-q/r],F)
def ocv(c):return oph(c.numerator())/oph(c.denominator())
def iev(f):return sum((CF(ocv(c))*l**int(ex[0])*k**(-int(ex[2])) for ex,c in f.dict().items() if int(ex[1])==0),CF.zero())
J=matrix(CF,6,6);J[0,3]=k/l;J[1,1]=1;J[2,0]=1;J[3,4]=-k;J[4,2]=1;J[5,5]=k*k
def c01(f):
 f=CF(f);nn=C(f.numerator());dd=C(f.denominator());assert dd[0],('kappa pole',f)
 return (F(nn[0]/dd[0]),F((nn[1]*dd[0]-nn[0]*dd[1])/dd[0]**2))
slots={tuple(map(int,e)):j for j,e in enumerate(two['normal_exponents'])};et2=(0,0,0,0,2,0);N0=[];N1=[];q0=[];q1=[]
for i,row in enumerate(inc['matrix'][:28].rows()):
 fac=k if i==0 or i>=10 else CF.one();nr=fac*vector(CF,[iev(f) for f in row])*J;cs=[c01(f) for f in nr];N0.append([x[0] for x in cs]);N1.append([x[1] for x in cs]);qq=c01(fac*iev(two['coefficient_matrix'][i,slots[et2]]));q0.append(qq[0]);q1.append(qq[1])
pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pa['ring'];Sx=S.gens();PK=S.base_ring();pph=PK.ring().hom([a,a*q/m,1/r,-q/r],F)
def pcv(c):return pph(c.numerator())/pph(c.denominator())
def ev(f):
 vals=[CF.one(),0,0,0,0,CF(l),l/k,0,1/k];return sum((CF(pcv(co))*prod(CF(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),CF.zero())
for it in pa['items']:
 f=it['full'];top=S({ex:co for ex,co in f.dict().items() if sum(ex)==4});lower=S({ex:co for ex,co in f.dict().items() if sum(ex)==3});mr=vector(CF,[ev(top.derivative(Sx[j])) for j in [1,2,3,4,6]]+[ev(lower)]);cs=[c01(k*f) for f in mr*J];N0.append([x[0] for x in cs]);N1.append([x[1] for x in cs]);qq=c01(k*ev(top.derivative(Sx[6]).derivative(Sx[6]))/2);q0.append(qq[0]);q1.append(qq[1])
N0=matrix(F,N0);N1=matrix(F,N1);assert N0==d['first_matrix'],'true transported first0';assert not N0.column(2) and not N1.column(2),'all-kappa E0 kernel first';rho,kap,E0=d['ring'].gens();oq0=vector(F,[f.monomial_coefficient(E0*E0) for f in d['rows']]);assert vector(F,q0)==oq0,'true eta² source coefficient0'
cols=[0,1,3,4,5];B0=N0.matrix_from_columns(cols);B1=N1.matrix_from_columns(cols);rr=d['range_rows'];A=d['range_matrix'];c2=vector(F,[f.monomial_coefficient(E0*E0) for f in d['range_compensation']]);v1=vector(F,q1)+B1*c2;c2k=-A.inverse()*vector(F,[v1[j] for j in rr]);resk=v1+B0*c2k;resrho=d['module_matrix'].column(0);M=matrix(F,[[resrho[j],resk[j]] for j in range(31)])
print('TRUE kappa range correction',list(c2k),flush=True);print('RANK rhoE0/kappaE0²',M.rank(),flush=True)
for j in range(31):
 if M[j,0] or M[j,1]:print('ROW',j,[str(c.factor()) if c else '0' for c in M.row(j)],flush=True)
if M.rank()==2:
 rp=M.transpose().pivots();print('FIRST TWO MINOR',rp,M.matrix_from_rows(rp).det().factor(),flush=True)
save({'field':F,'first_matrix':N0,'kappa_first_matrix':N1,'second_source_zero':q0,'second_source_kappa':q1,'range_kappa_quadratic':c2k,'kappa_quadratic_residual':resk,'rho_linear_residual':resrho,'module_matrix':M,'scope':'True complete original kappa*E0² source module, cyclic all material; whole ideal/amp wall pending'},str(p/'independent_subcritical_two_RR_inactive_pair_kappa_quadratic_module.sobj'),compress=False)
print('PASS true kappa mixed source',time.time()-st,flush=True)
