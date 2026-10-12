"""Complete original kappa=rho=0 paired boundary, generic base coefficients.

No normal truncation.  The three affine e1/e2/delta ranges are kept as
complete polynomials; boundary identities can be checked before elimination.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.time()
qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')
resource.setrlimit(resource.RLIMIT_AS,(vm+500*1024**2,vm+500*1024**2))
print('LinuxPID',os.getpid(),flush=True)
F=PolynomialRing(QQ,names=('a','r','H','Z','l')).fraction_field();a,r,H,Z,l=F.gens()
q=1+r;m=1+a*r;b=a*q/m;h=H/r;zz=q*Z/r
al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));sh=-2*H*(Z+1)/(Z*(H+1))
P=PolynomialRing(F,names=('s','v','E','e1','e2','delta'));s,v,E,e1,e2,de=P.gens()
vals=[P.one(),P.zero(),v,s,P.zero(),P(l),P(al*l),P.zero(),P.one()]
es=[E+sh*l,e1,e2];rows=[];labels=[];checks=[]
ph=qp['ring'].base_ring().ring().hom([a,r,H,Z],F);qc={}
def qcv(c):
 if c not in qc:qc[c]=ph(c.numerator())/ph(c.denominator())
 return qc[c]
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
QO=qp['ring'];off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset']
DC=matrix(QO,qp['D'])*matrix(QO,qp['C']);Dd=matrix(QO,qp['D'])*vector(QO,off)
rows.append(sum((pull(f.dict(),qcv)*es[j] for j,f in enumerate(qp['Q0_L_cofactors'])),P.zero())+pull(face(qp['full_Q0'],(2,1)),qcv));labels.append('Q0')
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for j,it in enumerate(ev['items']):
 rows.append(qcv(it['square_unit'])*es[int(it['square_component'])]**2+sum((pull(c.dict(),qcv)*es[k] for k,c in enumerate(it['A_p_quadratic'])),P.zero())+pull(it['pure_p4'].dict(),qcv));labels.append('even'+str(j))
for j,f in enumerate(qp['full_R_sources']):
 rows.append(sum((pull(DC[j,k].dict(),qcv)*es[k] for k in range(3)),P.zero())+pull(Dd[j].dict(),qcv)+de*pull(face(f,(1,2)),qcv));labels.append('R'+str(j))
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'));sp=load(str(p/'independent_supercritical_RR_small_support_sources.sobj'))
OK=next(iter(hi['items'][3]['H41_original'].values())).parent();oph=OK.ring().hom([a,b,h,zz],F);oc={}
def ocv(c):
 if c not in oc:oc[c]=oph(c.numerator())/oph(c.denominator())
 return oc[c]
for j,it in enumerate(hi['items']):
 rows.append(sum((pull(dd,ocv)*es[k] for k,dd in enumerate(it['L_cofactors_original'])),P.zero())+pull(it['H41_original'],ocv)+de*pull(sp['H13_original'][j],ocv));labels.append('H'+str(j))
 print('full H',j,'terms',len(rows[-1].dict()),'sec',round(time.time()-st,2),flush=True)
pair=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pair['ring'];Sx=S.gens();pk=S.base_ring();pph=pk.ring().hom([a,b,h,zz],F)
def pcv(c):return pph(c.numerator())/pph(c.denominator())
for j,it in enumerate(pair['items']):
 f22=S(face(it['full'],(2,2)))
 ff=(E+sh*l)*pull(f22.derivative(Sx[6]).dict(),pcv)+e1/(al*l)*pull(f22.derivative(Sx[1]).dict(),pcv)-e2/be*pull(f22.derivative(Sx[4]).dict(),pcv)+pull(face(it['full'],(3,1)),pcv)+de*pull(face(it['full'],(0,3)),pcv)
 rows.append(ff);labels.append('pair'+str(j))
 print('full pair',j,'terms',len(ff.dict()),flush=True)
zero={e1:0,e2:0,de:0};boundary=[P(f.subs(zero)) for f in rows]
for n,f in zip(labels,boundary):
 if f:print('NONZERO complete boundary',n,'terms',len(f.dict()),flush=True)
 else:checks.append('complete '+n+' vanishes for arbitrary s,v,E on e1=e2=delta=0')
save({'ring':P,'rows':rows,'labels':labels,'zero_range_boundary':boundary,'alpha':al,'beta':be,'shift':sh,'checks':checks,'scope':'Complete generic kappa=rho=0 paired source chart; no ideal or no-lift conclusion'},str(p/'independent_subcritical_two_RR_inactive_pair_generic_boundary.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'rows':len(rows),'seconds':time.time()-st};(p/'independent_subcritical_two_RR_inactive_pair_generic_boundary.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
