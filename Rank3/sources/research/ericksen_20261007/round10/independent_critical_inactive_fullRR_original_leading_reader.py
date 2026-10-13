"""Original-source classification of every critical inactive full-RR paired limit.
Use complete kappa=0 sources, delta a unit; F12 is not assumed separately zero.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(label,z):assert bool(z),label;checks.append(label)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));init=load(str(p/'independent_delta_subcritical_H_initials.sobj'));raw=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));blocks=load(str(p/'independent_inactive_RRtorus_normal_blocks.sobj'));bd=load(str(p/'independent_inactive_RRtorus_incoming_boundary.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
K=qf['ring'].base_ring();a,r,H,Z=K.gens();q=1+r;alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));shift=-2*H*(Z+1)/(Z*(H+1));P=PolynomialRing(K,names=('cl_l','cl_rho','cl_e0','cl_e1','cl_e2','cl_delta'));l,rho,e0,e1,e2,de=P.gens();es=(e0,e1,e2);imgs=(P(1),P(0),P(0),P(0),P(0),l,alpha*l,rho,P(1));F=P.fraction_field();oc=raw['coefficient_field'];mp=oc.ring().hom([a,a*q/(1+a*r),H/r,q*Z/r],K);cache={}
def cv(z):
 z=oc(z)
 if z not in cache:cache[z]=mp(z.numerator())/mp(z.denominator())
 return cache[z]
def pull(dd,old=False):return sum(((cv(co) if old else K(co))*prod(imgs[j]**int(ex[j]) for j in range(9)) for ex,co in dd.items() if not any(ex[j] and not imgs[j] for j in range(9))),P(0))
def part(f,bi):return {tuple(map(int,e)):co for e,co in f.dict().items() if (sum(e[:6]),sum(e[6:]))==bi}
DC=qf['D']*qf['C'];sources=[];labels=[]
Q0=sum((pull(f.dict())*es[j] for j,f in enumerate(qf['Q0_L_cofactors'])),P(0))+pull(part(qf['full_Q0'],(2,1)));sources.append(Q0);labels.append('Q0')
for j in (1,3):
 val=sum((pull(DC[j,k].dict())*es[k] for k in range(3)),P(0))+pull(part(qf['full_R_sources'][j],(4,0)))+de*pull(part(qf['full_R_sources'][j],(1,2)));sources.append(val);labels.append('R'+str(j))
pool=raw['polynomial_pool']
def face(rec,bi):
 out={}
 for sign,li,ri in rec['full_signed_product_DAG']:
  for ex,co in pool[li].items():
   if any(ex[k] for k in (1,2,3,4)):continue
   for ey,co2 in pool[ri].items():
    if any(ey[k] for k in (1,2,3,4)):continue
    ee=tuple(int(x)+int(y) for x,y in zip(ex,ey))
    if (sum(ee[:6]),sum(ee[6:]))==bi:out[ee]=out.get(ee,oc(0))+sign*co*co2
 return {e:co for e,co in out.items() if co}
for j in (1,13):
 val=sum((pull(cf,True)*es[k] for k,cf in enumerate(init['items'][j]['L_cofactors_original'])),P(0))+pull(face(raw['records'][j],(4,1)),True)+de*pull(face(raw['records'][j],(1,3)),True);sources.append(val);labels.append('H'+str(j));print('RAW H',j,'seconds',time.monotonic()-st,flush=True)
for i,f in enumerate(sources):
 ck('full critical affine source '+labels[i],all(sum(ex[2:])<=1 for ex in f.dict()))
 expectedrow=[6,1,3,8,20][i];h=bd['ring'].hom(P.gens(),P);ck('original source matches whole first boundary '+labels[i],f==h(bd['rows'][expectedrow]))
for i in range(1,5):ck('critical range/H e0 absent '+labels[i],sources[i].derivative(e0)==0)
M=matrix(P,[[f.derivative(v) for v in (e1,e2,de)] for f in sources[1:]]);seed=M.matrix_from_rows_and_columns([0,1],[0,1]).det();ck('protected two R e range nonzero',seed!=0)
pminor=M.matrix_from_rows([0,1,2]).det();qminor=M.matrix_from_rows([0,1,3]).det();PP=rho*rho/(beta*beta)-5*r*l*l-4*q;QQ=rho*rho/(beta*beta)+4*r*l*l+5*q;up=F(pminor)/(l**4*PP);uq=F(qminor)/(l**3*QQ)
ck('literal original P minor protected factor',up.numerator().total_degree()==0 and up.denominator().total_degree()==0);ck('literal original Q minor protected factor',uq.numerator().total_degree()==0 and uq.denominator().total_degree()==0);ck('exact canonical elimination difference',QQ-PP==9*(r*l*l+q));print('RANGE2 UNIT',K(F(seed/l**2)).factor(),flush=True);print('P UNIT',K(up).factor(),'Q UNIT',K(uq).factor(),flush=True)
# A delta-unit solution gives a nonzero vector (e1,e2,delta) in the kernel,
# hence both true minors vanish. Their protected units force PP=QQ=0.
# This proves l^2=-q/r,rho^2=-q beta^2 without using F12 separately.
U=PolynomialRing(K,names=('Lc','Rc'));ul,ur=U.gens();E=U.quotient(U.ideal(ul**2+q/r,ur**2+q*beta*beta),names=('L','R'));L,R=E.gens()
def redpol(f):return sum((co*L**int(ex[0])*R**int(ex[1])*0**sum(ex[2:]) for ex,co in P(f).dict().items()),E(0))
# Check the canonical delta coefficients vanish on all four roots, and that
# the resulting unique e graph is the correct generic shift, not cyclic e=0.
for i,f in enumerate(sources):
 cf=f.derivative(de);ck('all canonical critical delta coefficients zero '+labels[i],redpol(cf)==0)
 z=f.subs({e0:shift*l,e1:0,e2:0});z=z.subs({de:1});ck('canonical correct shifted graph satisfies '+labels[i],redpol(z)==0)
q0coef=Q0.derivative(e0);ck('Q0 third e pivot protected rho',F(q0coef/rho).numerator().total_degree()==0 and F(q0coef/rho).denominator().total_degree()==0 and q0coef!=0)
# Existing complete-source readers certify old base denominators. Here all
# newly used coefficient units are explicitly saved; none includes a material
# wall or a state amplitude beyond the chart l,rho units.
rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'true_critical_delta_unit_leading_classification':True,'no_separate_F12_zero_assumed':True,'all_four_root_e_graphs':True,'scope':__doc__}
save({'record':rec,'checks':checks,'ring':P,'sources':sources,'labels':labels,'matrix':M,'range_two_minor':seed,'P':PP,'Q':QQ,'P_unit':K(up),'Q_unit':K(uq),'Q0_e0':q0coef,'shift':shift},str(p/'independent_critical_inactive_fullRR_original_leading_reader.sobj'),compress=False);(p/'independent_critical_inactive_fullRR_original_leading_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,flush=True)
