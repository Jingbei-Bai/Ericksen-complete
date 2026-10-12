"""Original generic rho*E0 source and true range compensation."""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time();kd=load(str(p/'independent_subcritical_two_RR_inactive_pair_generic_kappa_module.sobj'));qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+400*1024**2,vm+400*1024**2))
print('LinuxPID',os.getpid(),flush=True)
d=kd['generic_boundary'];F=d['ring'].base_ring();a,r,H,Z,l=F.gens();q=1+r;m=1+a*r;P=PolynomialRing(F,'rho');rho=P.gen();vals=[P.one(),0,0,0,0,P(l),P(d['alpha']*l),rho,P.one()];checks=[]
def ck(n,t):
 if not t:raise AssertionError(n)
 checks.append(n)
ph=qp['ring'].base_ring().ring().hom([a,r,H,Z],F)
def qcv(c):return ph(c.numerator())/ph(c.denominator())
def pull(dd,cv):
 out=P.zero()
 for ee,co in dd.items():
  ex=tuple(map(int,ee));ex=ex+(0,)*(9-len(ex));mon=P.one()
  for j,n in enumerate(ex):
   if n:mon*=P(vals[j])**n
   if not mon:break
  if mon:out+=cv(co)*mon
 return out
QO=qp['ring'];DC=matrix(QO,qp['D'])*matrix(QO,qp['C']);raw=[pull(qp['Q0_L_cofactors'][0].dict(),qcv)]
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for it in ev['items']:raw.append(pull(it['A_p_quadratic'][0].dict(),qcv)+ (2*qcv(it['square_unit'])*d['shift']*l if int(it['square_component'])==0 else 0))
for j in range(6):raw.append(pull(DC[j,0].dict(),qcv))
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'));OK=next(iter(hi['items'][3]['H41_original'].values())).parent();oph=OK.ring().hom([a,a*q/m,H/r,q*Z/r],F)
def ocv(c):return oph(c.numerator())/oph(c.denominator())
for it in hi['items']:raw.append(pull(it['L_cofactors_original'][0],ocv))
pa=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pa['ring'];pk=S.base_ring();pph=pk.ring().hom([a,a*q/m,H/r,q*Z/r],F)
def pcv(c):return pph(c.numerator())/pph(c.denominator())
for it in pa['items']:
 f22=S({ex:co for ex,co in it['full'].dict().items() if (sum(ex[:6]),sum(ex[6:]))==(2,2)});raw.append(pull(f22.derivative(S.gen(6)).dict(),pcv))
for j,f in enumerate(raw):ck('true original E0 coefficient zero at rho0 '+str(j),not f[0])
qr=vector(F,[f[1] for f in raw]);N=kd['original_first_matrix'];B=N.matrix_from_columns([0,1,3,4,5]);rr=kd['range_rows'];A=kd['range_matrix'];cr=-A.inverse()*vector(F,[qr[j] for j in rr]);res=qr+B*cr
cyF=PolynomialRing(QQ,names=('a','r','l')).fraction_field();aa,rr0,ll=cyF.gens();cyph=F.ring().hom([aa,rr0,1,-1,ll],cyF)
def cyc(c):return cyph(c.numerator())/cyph(c.denominator())
old=load(str(p/'independent_subcritical_two_RR_inactive_pair_second_module.sobj'))
ON=old['first_matrix'].matrix_from_columns([0,1,3,4,5]);Oq=vector(cyF,[f.monomial_coefficient(old['ring'].gen(0)*old['ring'].gen(2)) for f in old['rows']]);Oa=ON.matrix_from_rows(rr);Oc=-Oa.inverse()*vector(cyF,[Oq[j] for j in rr]);Ores=Oq+ON*Oc
for j,c in enumerate(res):ck('independent cyclic true rhoE0 crossmatch '+str(j),cyc(c)==Ores[j])
M=matrix(F,[[res[j],kd['true_kappa_E0_squared_residual'][j],0] for j in range(31)]);fc=load(str(p/'independent_subcritical_two_RR_inactive_pair_boundary_fourth_certificate.sobj'));M[29,2]=fc['pair1_fourth'];minors=[]
for j in [5,11,13,23,25,30]:print('CYCLIC RHOE',j,cyc(res[j]).factor(),flush=True)
for i in [30]:
 for j in [17,19]:
  de=M.matrix_from_rows([i,j,29]).det();cv=cyc(de);ck('generic three-source generator determinant nonzero '+str((i,j)),bool(de));minors.append({'rows':[i,j,29],'determinant':de,'cyclic_value':cv});print('CYCLIC THREEGEN',i,j,cv.factor(),flush=True)
save({'generic_source':kd,'rho_source_polynomials':raw,'rho_second_range':cr,'true_rho_E0_residual':res,'three_generator_module':M,'three_generator_minors':minors,'checks':checks,'scope':'True generic original 3-generator module off paired amplitude/pivot walls; complete membership follows from boundary and whole-kernel identities'},str(p/'independent_subcritical_two_RR_inactive_pair_generic_rho_module.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_pair_generic_rho_module.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
