"""All-base/all-kappa original pure-axis first kernel, including full pairs."""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();kc=load(str(p/'independent_ratio_axis_normal_kernel_certificate.sobj'));inc=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+400*1024**2,vm+400*1024**2))
print('LinuxPID',os.getpid(),flush=True)
F=PolynomialRing(QQ,names=('a','r','H','Z','rho')).fraction_field();a,r,H,Z,rho=F.gens();q=1+r;m=1+a*r;be=H*(Z-1)/(Z*(H+1));C=PolynomialRing(F,'kap');k=C.gen();CF=C.fraction_field();OK=inc['coefficient_field'];ph=OK.ring().hom([a,a*q/m,H/r,q*Z/r],F);checks=[]
def cv(c):return ph(c.numerator())/ph(c.denominator())
def au(f):return sum((CF(cv(co))*(rho/k)**int(ex[0])*(1/k)**int(ex[1]) for ex,co in f.dict().items()),CF.zero())
bo=au(kc['beta']);ww=-k/be;ss=bo*ww;oldvec=vector(CF,[ss,0,0,ww,0,0])
def iev(f):return sum((CF(cv(co))*(rho/k)**int(ex[1])*(1/k)**int(ex[2]) for ex,co in f.dict().items() if int(ex[0])==0),CF.zero())
for j,row in enumerate(inc['matrix'].rows()):
 zz=vector(CF,[iev(f) for f in row])*oldvec;assert zz==0,('whole original first kernel',j);checks.append('whole generic all-kappa original first kernel '+str(j))
pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pa['ring'];PK=S.base_ring();pph=PK.ring().hom([a,a*q/m,H/r,q*Z/r],F)
def pcv(c):return pph(c.numerator())/pph(c.denominator())
def ev(f):
 vals=[CF.one(),0,0,0,0,0,0,rho/k,1/k];return sum((CF(pcv(co))*prod(CF(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),CF.zero())
for j,it in enumerate(pa['items']):
 top=S({ex:co for ex,co in it['full'].dict().items() if sum(ex)==4});zz=k*(ev(top.derivative(S.gen(1)))*ss+ev(top.derivative(S.gen(4)))*ww);assert not zz,('whole original full pair first kernel',j);checks.append('whole generic all-kappa FULL pair first kernel '+str(j))
L=PolynomialRing(F,names=('kap','l'));lk,ll=L.gens();LF=L.fraction_field();al=(Z-H)/(Z*(H+1));sh=-2*H*(Z+1)/(Z*(H+1))
def sev(f):
 vals=[LF.one(),0,0,0,0,LF(ll),ll*(al/lk+sh),rho/lk,1/lk];return sum((LF(pcv(co))*prod(LF(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),LF.zero())
for j,it in enumerate(pa['items']):
 top=S({ex:co for ex,co in it['full'].dict().items() if sum(ex)==4});assert sev(top)==0,('COMPLETE allbase allkappa pair S-prime',j);checks.append('COMPLETE allbase allkappa fullpair S-prime zero '+str(j))
cyF=PolynomialRing(QQ,names=('a','r','rho')).fraction_field();aa,rr,rrho=cyF.gens();cyph=F.ring().hom([aa,rr,1,-1,rrho],cyF);cyC=PolynomialRing(cyF,'kap');kk=cyC.gen();cyCF=cyC.fraction_field()
def cyc(f):
 f=CF(f)
 def cp(z):return cyC({int(e):cyph(c.numerator())/cyph(c.denominator()) for e,c in C(z).dict().items()})
 return cyCF(cp(f.numerator())/cp(f.denominator()))
assert cyc(ww)==-kk and cyc(ss)==(1+6*kk)/4;checks.append('cyclic regular first kernel s=(1+6kap)/4,w=-kap')
save({'coefficient_field':F,'kappa_ring':C,'kernel_s':ss,'kernel_w':ww,'old_normal_kernel':oldvec,'checks':checks,'scope':'Generic whole all-kappa pure-axis first kernel for all37 original sources AND all3 full pair sources'},str(p/'independent_subcritical_two_RR_inactive_axis_whole_kernel_reader.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_axis_whole_kernel_reader.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
