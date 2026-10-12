"""Exact entry into the paired chart, including ALL finite-e roots."""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();d=load(str(p/'independent_subcritical_two_RR_inactive_pair_generic_boundary.sobj'));qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+350*1024**2,vm+350*1024**2))
print('LinuxPID',os.getpid(),flush=True)
P=d['ring'];F=P.base_ring();a,r,H,Z,l=F.gens();s,v,E,e1,e2,de=P.gens();al=d['alpha'];be=d['beta'];sh=d['shift'];QO=qp['ring'];ph=QO.base_ring().ring().hom([a,r,H,Z],F);checks=[]
def ck(n,t):
 if not t:raise AssertionError(n)
 checks.append(n)
def cv(c):return ph(c.numerator())/ph(c.denominator())
C=PolynomialRing(F,names=('kap','rho','s','v','E','e1','e2'));k,rho,cs,cvv,cE,ce1,ce2=C.gens();CF=C.fraction_field();Rv=[al*l+k*(cE+sh*l),rho,CF.one()];pv=[CF.one(),(al/be*cs*rho+k*ce1)/Rv[0],cvv,cs,(rho*cvv-k*ce2)/be,CF(l)];vals=pv+Rv
def pull(f):return sum((CF(cv(co))*prod(CF(vals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),CF.zero())
LL=matrix(CF,[[pull(f) for f in row] for row in qp['L']]);ck('COMPLETE ORIGINAL L incidence exact, all normal values',LL*vector(CF,Rv)==k*vector(CF,[cE+sh*l,ce1,ce2]))
T=PolynomialRing(F,'RR1');tt=T.gen();leadvals=[T.one(),0,0,0,0,T(l),tt,0,T.one()]
def lp(f):return sum((T(cv(co))*prod(T(leadvals[j])**int(ex[j]) for j in range(9)) for ex,co in f.dict().items()),T.zero())
LLlead=matrix(T,[[lp(f) for f in row] for row in qp['L']]);ck('ORIGINAL leading L forces unique RR1=alpha*l',LLlead*vector(T,[tt,0,1])==vector(T,[tt-al*l,0,0]))
nv=[e1,e2,de];sv0={s:0,v:0};rf=[P(d['rows'][j].subs(sv0)) for j in [0,7,28]];A=matrix(F,[[F(f.derivative(x).constant_coefficient()) for x in nv] for f in rf]);
for i,f in enumerate(rf):ck('complete finite-e affine source '+str(i),f==sum((A[i,j]*nv[j] for j in range(3)),P.zero()))
ck('full finite-e affine determinant nonzero generically',bool(A.det()))
cyF=PolynomialRing(QQ,names=('a','r','l')).fraction_field();aa,rr,ll=cyF.gens();cyph=F.ring().hom([aa,rr,1,-1,ll],cyF)
def cyc(c):return cyph(c.numerator())/cyph(c.denominator())
ck('all-material/all-amplitude affine finite-e unit',cyc(A.det())==-9*(1+rr)**2*ll/(14336*rr**2))
fe=P(d['rows'][2].subs({s:0,v:0,e1:0,e2:0,de:0}));ue=F(fe.monomial_coefficient(E**2));ck('COMPLETE original even1 pure E², all finite roots',fe==ue*E**2)
ck('finite-e square coefficient explicit allbase unit',ue==2*H**2*Z**2/((Z-H)**2*(Z-1)**2));ck('finite-e cyclic unit',cyc(ue)==1/8)
pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pa['ring'];PK=S.base_ring();pph=PK.ring().hom([a,a*(1+r)/(1+a*r),H/r,(1+r)*Z/r],F)
zeroimgs=[CF(z.subs({k:0})) for z in vals]
for j,it in enumerate(pa['items']):
 zz=CF.zero()
 for ex,co in it['full'].dict().items():
  if (sum(ex[:6]),sum(ex[6:]))!=(2,2):continue
  cc=pph(co.numerator())/pph(co.denominator());zz+=CF(cc)*prod(zeroimgs[ii]**int(ex[ii]) for ii in range(9))
 ck('COMPLETE generic pair22 zero on full kappa0 incidence '+str(j),zz==0)
save({'generic_source':d,'L_source':qp['L'],'full_incidence_ring':C,'full_incidence_p':pv,'full_incidence_RR':Rv,'full_incidence_L':LL,'leading_L':LLlead,'finite_e_range_matrix':A,'finite_e_even1':fe,'finite_e_square_unit':ue,'checks':checks,'scope':'Every paired leading finite-e root lies on e=(shift*l,0,0), RRbar=(alpha*l,0,1), including both amplitude walls; full incidence chart exact and regular at lunit'},str(p/'independent_subcritical_two_RR_inactive_pair_finite_e_entry_certificate.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_pair_finite_e_entry_certificate.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
