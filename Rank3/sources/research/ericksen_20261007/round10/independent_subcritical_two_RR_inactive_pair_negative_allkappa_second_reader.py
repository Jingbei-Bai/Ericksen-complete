"""All-base/all-kappa negative-wall exact source second solution.
No matrix inverse, no frozen fiber: evaluate the COMPLETE original source
normal coefficient arrays on the candidate true second path, then reduce
every kappa coefficient in the protected quadratic material quotient.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();bd=load(str(p/'independent_subcritical_two_RR_inactive_pair_boundary_fourth_certificate.sobj'));inc=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));two=load(str(p/'independent_ratio_incoming12_generic_normal_two.sobj'));pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
F=bd['generic_boundary']['ring'].base_ring();a,r,H,Z,l=F.gens();q=1+r;m=1+a*r;C=PolynomialRing(F,'kap');k=C.gen();CF=C.fraction_field();OK=inc['coefficient_field'];oph=OK.ring().hom([a,a*q/m,H/r,q*Z/r],F);cc={};checks=[]
def ck(n,t):
 if not t:raise AssertionError(n)
 checks.append(n)
def ocv(c):
 if c not in cc:cc[c]=oph(c.numerator())/oph(c.denominator())
 return cc[c]
def iev(f):return sum((CF(ocv(co))*l**int(ex[0])*k**(-int(ex[2])) for ex,co in f.dict().items() if int(ex[1])==0),CF.zero())
slots={tuple(map(int,ex)):j for j,ex in enumerate(two['normal_exponents'])};eta2=(0,0,0,0,2,0);cs,cv=bd['true_second_range'];source=[]
for i,row in enumerate(inc['matrix'][:28].rows()):
 fac=k if i==0 or i>=10 else CF.one();f=fac*(iev(two['coefficient_matrix'][i,slots[eta2]])+iev(row[2])*cs+iev(row[1])*cv);source.append(CF(f));print('full kappa second source',i,'numterms',len(C(f.numerator()).dict()),'sec',round(time.time()-st,2),flush=True)
S=pa['ring'];sx=S.gens();PK=S.base_ring();pph=PK.ring().hom([a,a*q/m,H/r,q*Z/r],F);al=bd['generic_boundary']['alpha'];sh=bd['generic_boundary']['shift']
def pcv(c):return pph(c.numerator())/pph(c.denominator())
def ev(f):
 vals=[CF.one(),0,0,0,0,CF(l),l*(al/k+sh),0,1/k];return sum((CF(pcv(co))*prod(CF(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),CF.zero())
for j,it in enumerate(pa['items']):
 top=S({ex:co for ex,co in it['full'].dict().items() if sum(ex)==4});f=k*(ev(top.derivative(sx[6]).derivative(sx[6]))/2+ev(top.derivative(sx[3]))*cs+ev(top.derivative(sx[2]))*cv);source.append(CF(f))
B=PolynomialRing(QQ,names=('a','r','H','Z')).fraction_field();aa,rr,HH,ZZ=B.gens();U=PolynomialRing(B,'l');ll=U.gen();ph=F.ring().hom([aa,rr,HH,ZZ,ll],U);wall=ll**2+aa*rr/(1+aa*rr)
remainders=[]
for j,f in enumerate(source):
 nf=C(f.numerator());df=C(f.denominator());z=[]
 for ii,co in nf.dict().items():
  nr=ph(co.numerator())%wall;dr=ph(co.denominator())%wall;ck('allbase coefficient denominator unit at negative wall '+str((j,ii)),bool(dr));ck('COMPLETE all-kappa second numerator zero on negative wall '+str((j,ii)),nr==0);z.append((int(ii),nr,dr))
 remainders.append(z)
print('FULL original31 generic negative second zero',time.time()-st,flush=True)
save({'generic_boundary_source':bd,'allkappa_second_source_before_wall':source,'candidate_true_second_s_v':[cs,cv],'material_wall':wall,'exact_coefficient_remainders':remainders,'checks':checks,'scope':'COMPLETE original31 source true pure E0² solution on negative wall for ALL base/all kappa; no matrix rank specialization'},str(p/'independent_subcritical_two_RR_inactive_pair_negative_allkappa_second_reader.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'sources':31,'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_pair_negative_allkappa_second_reader.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
