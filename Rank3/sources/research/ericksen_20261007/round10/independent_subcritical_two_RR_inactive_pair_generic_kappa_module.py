"""Generic-base true kappa*E0^2 coefficients on the paired chart.

Complete original first/two coefficient archives are transported with the
actual shifted incoming family.  True boundary range correction is retained.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();bd=load(str(p/'independent_subcritical_two_RR_inactive_pair_boundary_fourth_certificate.sobj'));inc=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));two=load(str(p/'independent_ratio_incoming12_generic_normal_two.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
d=bd['generic_boundary'];P=d['ring'];F=P.base_ring();a,r,H,Z,l=F.gens();q=1+r;m=1+a*r;s,v,E,e1,e2,de=P.gens();al=d['alpha'];be=d['beta'];sh=d['shift'];C=PolynomialRing(F,'kap');k=C.gen();CF=C.fraction_field();checks=[]
def ck(n,t):
 if not t:raise AssertionError(n)
 checks.append(n)
OK=inc['coefficient_field'];oph=OK.ring().hom([a,a*q/m,H/r,q*Z/r],F);oc={}
def ocv(c):
 if c not in oc:oc[c]=oph(c.numerator())/oph(c.denominator())
 return oc[c]
def iev(f):return sum((CF(ocv(c))*l**int(ex[0])*k**(-int(ex[2])) for ex,c in f.dict().items() if int(ex[1])==0),CF.zero())
J=matrix(CF,6,6);J[0,3]=k/(l*(al+k*sh));J[1,1]=1;J[2,0]=1;J[3,4]=-k/be;J[4,2]=1;J[5,5]=k*k
def c01(f):
 f=CF(f);nn=C(f.numerator());dd=C(f.denominator());ck('regular kappa coefficient',bool(dd[0]))
 return (F(nn[0]/dd[0]),F((nn[1]*dd[0]-nn[0]*dd[1])/dd[0]**2))
slots={tuple(map(int,e)):j for j,e in enumerate(two['normal_exponents'])};et2=(0,0,0,0,2,0);N0=[];N1=[];q0=[];q1=[]
for i,row in enumerate(inc['matrix'][:28].rows()):
 fac=k if i==0 or i>=10 else CF.one();nr=fac*vector(CF,[iev(f) for f in row])*J;cs=[c01(f) for f in nr];N0.append([x[0] for x in cs]);N1.append([x[1] for x in cs]);qq=c01(fac*iev(two['coefficient_matrix'][i,slots[et2]]));q0.append(qq[0]);q1.append(qq[1])
 ck('whole all-kappa original first E0 kernel '+str(i),not nr[2])
pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pa['ring'];Sx=S.gens();PK=S.base_ring();pph=PK.ring().hom([a,a*q/m,H/r,q*Z/r],F)
def pcv(c):return pph(c.numerator())/pph(c.denominator())
def ev(f):
 vals=[CF.one(),0,0,0,0,CF(l),l*(al/k+sh),0,1/k];return sum((CF(pcv(co))*prod(CF(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),CF.zero())
for j,it in enumerate(pa['items']):
 f=it['full'];top=S({ex:co for ex,co in f.dict().items() if sum(ex)==4});lower=S({ex:co for ex,co in f.dict().items() if sum(ex)==3});mr=vector(CF,[ev(top.derivative(Sx[z])) for z in [1,2,3,4,6]]+[ev(lower)]);nr=k*mr*J;cs=[c01(f) for f in nr];N0.append([x[0] for x in cs]);N1.append([x[1] for x in cs]);qq=c01(k*ev(top.derivative(Sx[6]).derivative(Sx[6]))/2);q0.append(qq[0]);q1.append(qq[1]);ck('whole all-kappa pair first E0 kernel '+str(j),not nr[2])
N0=matrix(F,N0);N1=matrix(F,N1);orig=matrix(F,[[F(f.derivative(x).constant_coefficient()) for x in [s,v,E,e1,e2,de]] for f in d['rows']]);ck('generic original normalized first boundary crossmatch',N0==orig)
c2=vector(F,list(bd['true_second_range'])+[0,0,0]);cols=[0,1,3,4,5];B0=N0.matrix_from_columns(cols);B1=N1.matrix_from_columns(cols);rr=[0,1,2,7,28];A=B0.matrix_from_rows(rr);v1=vector(F,q1)+B1*c2;c2k=-A.inverse()*vector(F,[v1[j] for j in rr]);resk=v1+B0*c2k
print('generic true kappa range correction nonzero',[j for j,c in enumerate(c2k) if c],flush=True)
for j in [11,13,17,19,23,25,28,29,30]:print('kappa residual',j,'nonzero',bool(resk[j]),flush=True)
cyF=PolynomialRing(QQ,names=('a','r','l')).fraction_field();aa,rr0,ll=cyF.gens();cyph=F.ring().hom([aa,rr0,1,-1,ll],cyF)
def cyc(c):return cyph(c.numerator())/cyph(c.denominator())
old=load(str(p/'independent_subcritical_two_RR_inactive_pair_kappa_quadratic_module.sobj'))
for j,c in enumerate(resk):ck('independent cyclic generic-kappa crossmatch '+str(j),cyc(c)==old['kappa_quadratic_residual'][j])
save({'generic_boundary':d,'original_first_matrix':N0,'kappa_first_matrix':N1,'source_eta2_zero':q0,'source_eta2_kappa':q1,'range_rows':rr,'range_matrix':A,'true_kappa_second_range':c2k,'true_kappa_E0_squared_residual':resk,'checks':checks,'scope':'Generic complete source kappa E0² residual and whole all-kappa E0 first kernel; whole legal saturation proof pending'},str(p/'independent_subcritical_two_RR_inactive_pair_generic_kappa_module.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_pair_generic_kappa_module.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
