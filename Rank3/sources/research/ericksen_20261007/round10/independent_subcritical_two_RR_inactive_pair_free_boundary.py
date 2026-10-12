"""Complete cyclic kappa=rho=delta=e1=e2=0 paired source scheme."""
from sage.all import *
from pathlib import Path
import os,resource,time
p=Path(__file__).resolve().parent;st=time.time();qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
F=PolynomialRing(QQ,names=('a','r','l')).fraction_field();a,r,l=F.gens();q=1+r;m=1+a*r;P=PolynomialRing(F,names=('s','v','E'));s,v,E=P.gens();vals=[P.one(),P.zero(),v,s,P.zero(),P(l),P(l),P.zero(),P.one()]
ph=qp['ring'].base_ring().ring().hom([a,r,1,-1],F)
def qcv(c):return ph(c.numerator())/ph(c.denominator())
def pull(dd,cv):
 out=P.zero()
 for ee,co in dd.items():
  ex=tuple(map(int,ee));ex=ex+(0,)*(9-len(ex));mon=P.one()
  for j,n in enumerate(ex):
   if n:mon*=vals[j]**n
   if not mon:break
  if mon:out+=cv(co)*mon
 return out
def face(f,ij):return {tuple(map(int,ee)):co for ee,co in f.dict().items() if (sum(ee[:6]),sum(ee[6:]))==ij}
QO=qp['ring'];off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset'];DC=matrix(QO,qp['D'])*matrix(QO,qp['C']);Dd=matrix(QO,qp['D'])*vector(QO,off);rows=[];labels=[]
rows.append(pull(qp['Q0_L_cofactors'][0].dict(),qcv)*E+pull(face(qp['full_Q0'],(2,1)),qcv));labels.append('Q0')
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for j,it in enumerate(ev['items']):
 rows.append((qcv(it['square_unit'])*E*E if int(it['square_component'])==0 else P.zero())+pull(it['A_p_quadratic'][0].dict(),qcv)*E+pull(it['pure_p4'].dict(),qcv));labels.append('even'+str(j))
for j in range(6):rows.append(pull(DC[j,0].dict(),qcv)*E+pull(Dd[j].dict(),qcv));labels.append('R'+str(j))
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'));OK=next(iter(hi['items'][3]['H41_original'].values())).parent();oph=OK.ring().hom([a,a*q/m,1/r,-q/r],F)
def ocv(c):return oph(c.numerator())/oph(c.denominator())
for j,it in enumerate(hi['items']):rows.append(pull(it['L_cofactors_original'][0],ocv)*E+pull(it['H41_original'],ocv));labels.append('H'+str(j))
pair=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pair['ring'];Sx=S.gens();pk=S.base_ring();pkph=pk.ring().hom([a,a*q/m,1/r,-q/r],F)
def pcv(c):return pkph(c.numerator())/pkph(c.denominator())
for j,it in enumerate(pair['items']):f22=S(face(it['full'],(2,2)));rows.append(E*pull(f22.derivative(Sx[6]).dict(),pcv)+pull(face(it['full'],(3,1)),pcv));labels.append('pair'+str(j))
for n,f in zip(labels,rows):
 if f:print('FULL_BOUNDARY',n,f.factor(),flush=True)
save({'ring':P,'rows':rows,'labels':labels,'scope':'Complete original paired cyclic boundary source scheme; free sheet membership pending'},str(p/'independent_subcritical_two_RR_inactive_pair_free_boundary.sobj'),compress=False)
print('SAVED',len(rows),time.time()-st,flush=True)
