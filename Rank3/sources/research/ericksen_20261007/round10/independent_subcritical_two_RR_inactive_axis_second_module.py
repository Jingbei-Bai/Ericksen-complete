from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.time();F=PolynomialRing(QQ,names=('ax_a','ax_r','ax_rho')).fraction_field();a,r,rho=F.gens();q=1+r;m=1+a*r
P=PolynomialRing(F,names=('l','kap','s','z','w','E0','e1','delta'));l,kap,s,z,w,E0,e1,de=P.gens();norms=[w,E0,e1,de];pars=[l,kap,s,z];checks=[]
def ck(n,b):
 if not b:raise AssertionError(n)
 checks.append(n)
def tr(f,d=2):return P({ee:co for ee,co in P(f).dict().items() if sum(ee)<=d})
pv=[P.one(),s,w/rho+kap*z/rho,l*s/rho-kap*e1/rho+kap*E0*s/rho,w,l];Rv=[l+kap*E0,P(rho),P.one()];vals=pv+Rv;es=[E0,e1,z]
qp=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));QO=qp['ring'];QK=QO.base_ring();qph=QK.ring().hom([a,r,1,-1],F);qcache={}
def qcv(c):
 if c not in qcache:qcache[c]=qph(c.numerator())/qph(c.denominator())
 return qcache[c]
old=load(str(p/'independent_delta_H_global_L_module.sobj'))['original_field'];oph=old.ring().hom([a,a*q/m,1/r,-q/r],F);ocache={}
def ocv(c):
 if c not in ocache:ocache[c]=oph(c.numerator())/oph(c.denominator())
 return ocache[c]
def pull(dd,cv,d=2):
 out=P.zero()
 for ee,co in dd.items():
  e=tuple(map(int,ee));e=e+(0,)*(9-len(e));mindeg=e[1]+e[2]+2*e[3]+e[4]+e[5]+e[6]
  if mindeg>d:continue
  mon=P.one()
  for j,n in enumerate(e):
   if n:mon=tr(mon*tr(vals[j]**n,d),d)
   if not mon:break
  if mon:out+=P(cv(co))*mon
 return tr(out,d)
def face(f,ij):return {tuple(map(int,ee)):co for ee,co in f.dict().items() if (sum(ee[:6]),sum(ee[6:]))==ij}
DC=matrix(QO,qp['D'])*matrix(QO,qp['C']);off=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'))['offset'];Dd=matrix(QO,qp['D'])*vector(QO,off)
rows=[];labels=[]
Q=qp['full_Q0'];q0=sum((tr(pull(f.dict(),qcv)*es[j]) for j,f in enumerate(qp['Q0_L_cofactors'])),P.zero())+pull(face(Q,(2,1)),qcv)+tr(kap*pull(face(Q,(3,0)),qcv))+tr(kap*de*pull(face(Q,(0,2)),qcv));rows.append(tr(q0));labels.append(('Q0',0))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
for j,it in enumerate(ev['items']):
 ff=P(qcv(it['square_unit']))*es[int(it['square_component'])]**2+sum((tr(pull(c.dict(),qcv)*es[h]) for h,c in enumerate(it['A_p_quadratic'])),P.zero())+pull(it['pure_p4'].dict(),qcv)+tr(kap*de*pull(face(it['full_source'],(2,1)),qcv));rows.append(tr(ff));labels.append(('even',j))
for j,f in enumerate(qp['full_R_sources']):
 ff=sum((tr(pull(DC[j,h].dict(),qcv)*es[h]) for h in range(3)),P.zero())+pull(Dd[j].dict(),qcv)+tr(de*pull(face(f,(1,2)),qcv))+tr(kap*de*pull(face(f,(2,1)),qcv));rows.append(tr(ff));labels.append(('R',j))
hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'));sp=load(str(p/'independent_supercritical_RR_small_support_sources.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));pool=full['polynomial_pool']
vm=int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE');resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
# Only a one-normal H50 and the zero-normal H22 value enter the second jet.
def raw_small_face(idx,bi,deg):
 out={}
 for sg,li,ri in full['records'][idx]['full_signed_product_DAG']:
  for ex,c in pool[li].items():
   bx=(sum(ex[:6]),sum(ex[6:]));mx=sum(ex[h] for h in [1,2,3,4,5,6])
   if mx>deg or bx[0]>bi[0] or bx[1]>bi[1]:continue
   for ey,d in pool[ri].items():
    if sum(ey[h] for h in [1,2,3,4,5,6])+mx>deg:continue
    if (bx[0]+sum(ey[:6]),bx[1]+sum(ey[6:]))!=bi:continue
    ee=tuple(int(x)+int(y) for x,y in zip(ex,ey));out[ee]=out.get(ee,old.zero())+sg*c*d
 return {ee:co for ee,co in out.items() if co}
for j,it in enumerate(hi['items']):
 ff=sum((tr(pull(di,ocv)*es[h]) for h,di in enumerate(it['L_cofactors_original'])),P.zero())+pull(it['H41_original'],ocv)+tr(de*pull(sp['H13_original'][j],ocv))+tr(kap*pull(raw_small_face(j,(5,0),1),ocv))+tr(kap*de*pull(raw_small_face(j,(2,2),0),ocv));rows.append(tr(ff));labels.append(('H',j));print('jet H',j,'seconds',round(time.time()-st,2),flush=True)
pair=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'));pairk=pair['ring'].base_ring();pkph=pairk.ring().hom([a,a*q/m,1/r,-q/r],F)
def pcv(c):return pkph(c.numerator())/pkph(c.denominator())
for j,it in enumerate(pair['items']):
 lead=pull(face(it['full'],(2,2)),pcv,3);ck('pair22 true division '+str(j),all(int(ee[1])>=1 for ee in lead.dict()));ff=tr(lead/kap)+pull(face(it['full'],(3,1)),pcv)+tr(kap*pull(face(it['full'],(4,0)),pcv))+tr(de*pull(face(it['full'],(0,3)),pcv))+tr(kap*de*pull(face(it['full'],(1,2)),pcv));rows.append(tr(ff));labels.append(('pair',j))
for j,f in enumerate(rows):
 ck('whole Sline jet '+str(j),f.subs({s:0,z:0,w:0,E0:0,e1:0,de:0})==0)
 ck('first kernel sz '+str(j),tr(f,1).derivative(s)==tr(f,1).derivative(z)==0)
N=matrix(F,[[f.derivative(v).constant_coefficient() for v in norms] for f in rows]);ck('first range rank4',N.rank()==4)
probe=load(str(p/'independent_subcritical_two_RR_inactive_incoming_normal_probe.sobj'));Ap=probe['axis_matrix'];pk=Ap.base_ring().base_ring() if False else None
# Two protected H range charts, corresponding to the already exhaustive
# original critical R3/Q0/H9 versus R3/Q0/H15 delta pivots.
items=[]
generators=[l*s,l*z,kap*s,kap*z,s*s,s*z,z*z]
for rr in [(0,1,7,19),(0,1,7,25)]:
 A=N.matrix_from_rows(rr);adet=A.det();print('RANGE',rr,adet.factor() if adet else 0,flush=True)
 if not adet:continue
 qrows=[tr(f.subs({v:0 for v in norms}),2) for f in rows];comp=-A.inverse()*vector(P,[qrows[j] for j in rr]);res=[tr(qrows[j]+sum((N[j,h]*comp[h] for h in range(4)),P.zero()),2) for j in range(len(rows))]
 for j,f in enumerate(res):ck('residual seven-generator span '+str(rr)+'/'+str(j),f==sum((f.monomial_coefficient(g)*g for g in generators),P.zero()))
 MM=matrix(F,[[f.monomial_coefficient(g) for g in generators] for f in res]);rank=MM.rank();piv=MM.transpose().pivots();print('MODULE',rr,'rank',rank,'pivots',piv,flush=True)
 if rank==7:print('MINOR7',rr,MM.matrix_from_rows(piv).det().factor(),flush=True)
 else:
  print('KERNEL',rr,[list(v) for v in MM.right_kernel().basis()],flush=True)
 items.append({'range_rows':rr,'range_matrix':A,'range_det':adet,'range_compensation':comp,'residuals':res,'module_matrix':MM,'module_rank':rank,'module_pivots':piv})
save({'ring':P,'rows':rows,'labels':labels,'first_range_matrix':N,'normal_order':norms,'kernel_order':[s,z],'generators':generators,'items':items,'checks':checks,'scope':'Exact cyclic all-material source second module diagnostic; full-relative membership and generic-base audit required before whole theorem'},str(p/'independent_subcritical_two_RR_inactive_axis_second_module.sobj'),compress=False)
rep={'PASS':True,'checks':len(checks),'seconds':time.time()-st,'module_ranks':[it['module_rank'] for it in items]};(p/'independent_subcritical_two_RR_inactive_axis_second_module.json').write_text(json.dumps(rep,indent=2));print(rep,flush=True)
