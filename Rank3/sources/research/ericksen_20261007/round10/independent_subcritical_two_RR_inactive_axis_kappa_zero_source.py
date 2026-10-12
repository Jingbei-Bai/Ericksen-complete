"""Complete kappa=0 source chart for a legitimate kappa-saturation test."""
from sage.all import *
from pathlib import Path
import os,resource,time,sys
p=Path(__file__).resolve().parent;st=time.time();mode=sys.argv[1] if len(sys.argv)>1 else 'cyclic'
qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
print('LinuxPID',os.getpid(),'mode',mode,flush=True)
if mode=='adapted':
 F=PolynomialRing(QQ,names=('a','r','H','Z','rho')).fraction_field();a,r,H,Z,rho=F.gens();q=1+r;m=1+a*r;b=a*q/m;h=H/r;zz=q*Z/r
elif mode=='generic':
 F=PolynomialRing(QQ,names=('a','b','h','z','rho')).fraction_field();a,b,h,zz,rho=F.gens();r=(a-b)/(a*(b-1));q=1+r;m=1+a*r;H=r*h;Z=r*zz/q
else:
 F=PolynomialRing(QQ,names=('a','r','rho')).fraction_field();a,r,rho=F.gens();q=1+r;m=1+a*r;b=a*q/m;H=F.one();Z=-F.one();h=H/r;zz=q*Z/r
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));shift=-2*H*(Z+1)/(Z*(H+1))
P=PolynomialRing(F,names=('l','s','w','E0','e1','e2','delta'));l,s,w,E0,e1,e2,de=P.gens();pv=[P.one(),s,beta*w/rho,beta*l*s/rho,w,l];Rv=[alpha*l,P(rho),P.one()];vals=pv+Rv;es=[E0+shift*l,e1,e2]
QO=qp['ring'];QK=QO.base_ring();qph=QK.ring().hom([a,r,H,Z],F);qc={}
def qcv(c):
 if c not in qc:qc[c]=qph(c.numerator())/qph(c.denominator())
 return qc[c]
def pull(dd,cv):
 out=P.zero()
 for ee,c in dd.items():
  ex=tuple(map(int,ee));ex=ex+(0,)*(9-len(ex));mon=P.one()
  for j,n in enumerate(ex):
   if n:mon*=vals[j]**n
   if not mon:break
  if mon:out+=cv(c)*mon
 return out
def face(f,ij):return {tuple(map(int,ee)):co for ee,co in f.dict().items() if (sum(ee[:6]),sum(ee[6:]))==ij}
off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset'];DC=matrix(QO,qp['D'])*matrix(QO,qp['C']);Dd=matrix(QO,qp['D'])*vector(QO,off)
rows={};rows['Q0']=sum((pull(f.dict(),qcv)*es[j] for j,f in enumerate(qp['Q0_L_cofactors'])),P.zero())+pull(face(qp['full_Q0'],(2,1)),qcv)
for j in ([1,3] if mode in ['generic','adapted'] else range(6)):rows['R'+str(j)]=sum((pull(DC[j,h].dict(),qcv)*es[h] for h in range(3)),P.zero())+pull(Dd[j].dict(),qcv)+de*pull(face(qp['full_R_sources'][j],(1,2)),qcv)
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for j,it in enumerate(ev['items']):rows['even'+str(j)]=P(qcv(it['square_unit']))*es[int(it['square_component'])]**2+sum((pull(c.dict(),qcv)*es[h] for h,c in enumerate(it['A_p_quadratic'])),P.zero())+pull(it['pure_p4'].dict(),qcv)
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'));sp=load(str(p/'independent_supercritical_RR_small_support_sources.sobj'));ok=sp['coefficient_field'] if 'coefficient_field' in sp else None
# Original H coefficient field is recovered from the genuine dictionaries.
cc=next(iter(hi['items'][3]['H41_original'].values()));OK=cc.parent();oph=OK.ring().hom([a,b,h,zz],F);oc={}
def ocv(c):
 if c not in oc:oc[c]=oph(c.numerator())/oph(c.denominator())
 return oc[c]
for j in ([3,7,9,15] if mode in ['generic','adapted'] else range(18)):
 it=hi['items'][j];rows['H'+str(j)]=sum((pull(dd,ocv)*es[h] for h,dd in enumerate(it['L_cofactors_original'])),P.zero())+pull(it['H41_original'],ocv)+de*pull(sp['H13_original'][j],ocv)
 print('full H',j,'terms',len(rows['H'+str(j)].dict()),'sec',round(time.time()-st,2),flush=True)
pair=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));S=pair['ring'];Sx=S.gens();pk=S.base_ring();pkph=pk.ring().hom([a,b,h,zz],F)
def pcv(c):return pkph(c.numerator())/pkph(c.denominator())
for j,it in (enumerate(pair['items']) if mode not in ['generic','adapted'] else []):
 f22=S(face(it['full'],(2,2)))
 ff=(E0+shift*l)*pull(f22.derivative(Sx[6]).dict(),pcv)+e2/rho*pull(f22.derivative(Sx[2]).dict(),pcv)+beta/(alpha*rho)*((E0+shift*l)*s-e1)*pull(f22.derivative(Sx[3]).dict(),pcv)+pull(face(it['full'],(3,1)),pcv)+de*pull(face(it['full'],(0,3)),pcv)
 rows['pair'+str(j)]=ff;print('full pair',j,'terms',len(ff.dict()),flush=True)
save({'ring':P,'rows':rows,'mode':mode,'alpha':alpha,'beta':beta,'shift':shift},str(p/('independent_subcritical_two_RR_inactive_axis_kappa_zero_'+mode+'.sobj')),compress=False)
print('SAVED complete kappa0 source',len(rows),'sec',time.time()-st,flush=True)
if mode=='adapted':sys.exit(0)
# Three true affine range equations in E0,e1,delta; all denominator guards retained.
PF=P.fraction_field();names=['Q0','R3','H9'];nv=[E0,e1,de];A=matrix(PF,[[rows[n].derivative(v) for v in nv] for n in names]);rhs=vector(PF,[rows[n].subs({v:0 for v in nv}) for n in names])
for n in names:assert all(rows[n].derivative(v).derivative(vv)==0 for v in nv for vv in nv)
sol=-A.inverse()*rhs;print('AFFINE RANGE solved',time.time()-st,flush=True);sub=dict(zip(nv,sol));reds={n:PF(f.subs(sub)) for n,f in rows.items()};print('H3 reduced-zero',not reds['H3'],'terms',len(reds['H3'].numerator().dict()),flush=True)
out={'ring':P,'rows':rows,'range_matrix':A,'range_solution':sol,'reduced':reds,'mode':mode}
save(out,str(p/('independent_subcritical_two_RR_inactive_axis_kappa_zero_'+mode+'_reduced.sobj')),compress=False)
if mode=='cyclic':
 print('H3 numerator factor',reds['H3'].numerator().factor() if reds['H3'] else 0,flush=True)
 print('H3 denominator factor',reds['H3'].denominator().factor(),flush=True)
 print('even numer terms',[len(reds['even'+str(j)].numerator().dict()) for j in range(3)],flush=True)
print('PASS exact complete kappa0 reduction',time.time()-st,flush=True)
