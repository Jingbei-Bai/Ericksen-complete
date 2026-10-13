"""Direct general-base finite-beta original-source sigma/quadratic pullbacks.
Every original H source and complete even RR-linear module is retained.
Four canonical root embeddings are kept; no frozen-beta extension is assumed.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(label,z):assert bool(z),label;checks.append(label)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));hi=load(str(p/'independent_delta_subcritical_H_initials.sobj'));hb=load(str(p/'independent_supercritical_first_boundary_probe.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));em=load(str(p/'independent_delta_even_L_module.sobj'));kv=load(str(p/'independent_supercritical_inactive_generic_resonance_kernel.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+850*1024**2,)*2)
K=PolynomialRing(QQ,names=('bd_r','bd_H','bd_Z','bd_beta')).fraction_field();r,H,Z,beta=K.gens();q=1+r;a=-q/r;m=-r;alpha=(Z-H)/(Z*(H+1));Bbase=H*(Z-1)/(Z*(H+1));d=-2*H*(Z+1)/(Z*(H+1));lc=-q/r;bc=-q
U=PolynomialRing(K,names=('bd_C','bd_B'));uc,ub=U.gens();E=U.quotient(U.ideal(uc**2-lc,ub**2-bc),names=('C','B'));c,b=E.gens();bases=(E(1),c,b,c*b)
def ec(f):return vector(K,[f.lift().monomial_coefficient(v) for v in (U(1),uc,ub,uc*ub)])
inv_cache={}
def iv(f):
 f=E(f)
 if f not in inv_cache:
  A=matrix(K,4,4,lambda i,j:ec(f*bases[j])[i]);v=A.solve_right(vector(K,[1,0,0,0]));z=sum((v[j]*bases[j] for j in range(4)),E(0));ck('root inverse '+str(len(inv_cache)),f*z==1);inv_cache[f]=z
 return inv_cache[f]
Q=PolynomialRing(K,names=('bd_C','bd_B','bd_t'));xc,xb,t0=Q.gens();T=Q.quotient(Q.ideal(xc**2-lc,xb**2-bc,t0**3),names=('C','B','T'));ct,bt,t=T.gens();uh=U.hom([ct,bt],T)
def eb(f):return uh(E(f).lift())
def coeff_t(f,n):return sum((co*c**int(ex[0])*b**int(ex[1]) for ex,co in T(f).lift().dict().items() if int(ex[2])==n),E(0))
qK=qf['ring'].base_ring();qp=qK.ring().hom([a,r,H,Z],K);old=full['coefficient_field'];op=old.ring().hom([a,a*q/m,H/r,q*Z/r],K);qc={};oc={}
def cv(co,original=False):
 co=(old if original else qK)(co);cache=oc if original else qc;ph=op if original else qp
 if co not in cache:cache[co]=ph(co.numerator())/ph(co.denominator())
 return cache[co]
def pull(dd,imgs,original=False,ring=T):
 out=ring(0)
 for ex,co in dd.items():
  vv=ring(1)
  for j,nn in enumerate(ex):
   if nn:vv*=imgs[j]**int(nn)
   if not vv:break
  if vv:out+=cv(co,original)*vv
 return out
def part(f,bi):return {tuple(map(int,e)):co for e,co in f.dict().items() if (sum(e[:6]),sum(e[6:]))==bi}
def der(dd,j):
 out={}
 for ex,co in dd.items():
  ee=list(map(int,ex))
  if ee[j]:co*=ee[j];ee[j]-=1;out[tuple(ee)]=co
 return out
DC=qf['D']*qf['C'];pbase=(E(1),E(0),E(0),E(0),E(0),c,alpha*c,Bbase*b,E(1));F12=[part(f,(1,2)) for f in qf['full_R_sources']]
def euler(dd,inds):return {ex:co*sum(ex[j] for j in inds) for ex,co in dd.items() if sum(ex[j] for j in inds)}
A=matrix(E,[[pull(euler(F12[i],inds),pbase,ring=E) for inds in ((1,4,7),(3,5,6))] for i in (1,3)]);Ai=matrix(E,[[A[1,1],-A[0,1]],[-A[1,0],A[0,0]]])*iv(A.det());ck('true original two R range inverse',A*Ai==identity_matrix(E,2));SR=vector(E,[pull(part(qf['full_R_sources'][i],(2,1)),pbase,ring=E)+c*d*pull(der(F12[i],6),pbase,ring=E) for i in (1,3)]);siglog=-Ai*SR
vph=kv['field'].ring().hom([r,H,Z],K);nu=[vph(z.numerator())/vph(z.denominator()) for z in kv['kernel']];n=[z*t for z in nu]
def path(bb,cc,ns,sigma):
 s,w,e0,e1,e2=ns;kappa=sigma*beta;es=(cc*(d+e0),bb*cc*e1,bb*e2)
 deninv=1/alpha*(1-kappa*(d+e0)/alpha+kappa**2*(d+e0)**2/alpha**2)
 imgs=(T(1),bb*(alpha*s+kappa*e1)*deninv,w+kappa*e2/Bbase,cc*s,bb*w,cc,cc*(alpha+kappa*(d+e0)),Bbase*bb,T(1))
 return imgs,es
initial,es0=path(bt,ct,n,T(0))
def GR(i,imgs,es):return sum((pull(DC[i,j].dict(),imgs)*es[j] for j in range(3)),T(0))+pull(part(qf['full_R_sources'][i],(4,0)),imgs)
GRquad=vector(E,[coeff_t(GR(i,initial,es0),2) for i in (1,3)]);ck('kernel true R first beta coefficient zero',all(coeff_t(GR(i,initial,es0),1)==0 for i in (1,3)));quadlog=-Ai*GRquad
sigimgs,siges=path(bt*(1+beta*eb(siglog[0])*t),ct*(1+beta*eb(siglog[1])*t),[T(0)]*5,t)
quadimgs,quades=path(bt*(1+beta*eb(quadlog[0])*t*t),ct*(1+beta*eb(quadlog[1])*t*t),n,T(0))
for j,i in enumerate((1,3)):
 sigR=pull(F12[i],sigimgs)+beta*GR(i,sigimgs,siges)+t*beta*pull(part(qf['full_R_sources'][i],(2,1)),sigimgs);ck('true original R sigma range path '+str(i),coeff_t(sigR,0)==0 and coeff_t(sigR,1)==0)
 quadR=pull(F12[i],quadimgs)+beta*GR(i,quadimgs,quades);ck('true original R quadratic range path '+str(i),quadR==0)
print('R RANGE PATHS READY',time.monotonic()-st,flush=True)
sigma=[];quadratic=[];labels=[]
def add_source(label,fs,fq):
 ck('whole-beta source value sigma '+label,coeff_t(fs,0)==0);ck('whole-beta kernel value+linear '+label,coeff_t(fq,0)==0 and coeff_t(fq,1)==0)
 z=coeff_t(fs,1);w=coeff_t(fq,2)
 for j,v in enumerate(ec(z)):ck('all finite beta sigma coefficient '+label+':'+str(j),v.derivative(beta)==0)
 for j,v in enumerate(ec(w)):ck('all finite beta kernel quadratic '+label+':'+str(j),v.derivative(beta)==0)
 labels.append(label);sigma.append(z);quadratic.append(w)
def GQ(imgs,es):return sum((pull(ff.dict(),imgs)*es[j] for j,ff in enumerate(qf['Q0_L_cofactors'])),T(0))+pull(part(qf['full_Q0'],(2,1)),imgs)
add_source('Q0',GQ(sigimgs,siges)+t*pull(part(qf['full_Q0'],(0,2)),sigimgs)+t*beta*pull(part(qf['full_Q0'],(3,0)),sigimgs),GQ(quadimgs,quades))
pool=full['polynomial_pool'];h22cache={}
def face_pair(rec,bi):
 out={}
 for sign,li,ri in rec['full_signed_product_DAG']:
  for ex,aa0 in pool[li].items():
   if any(ex[j] for j in (1,2,3,4)):continue
   for ey,bb0 in pool[ri].items():
    if any(ey[j] for j in (1,2,3,4)):continue
    ee=tuple(int(x)+int(y) for x,y in zip(ex,ey))
    if (sum(ee[:6]),sum(ee[6:]))==bi:out[ee]=out.get(ee,old(0))+sign*aa0*bb0
 return {ee:co for ee,co in out.items() if co}
for i,it in enumerate(hi['items']):
 def GH(imgs,es):return sum((pull(ff,imgs,True)*es[j] for j,ff in enumerate(it['L_cofactors_original'])),T(0))+pull(it['H41_original'],imgs,True)
 h13=hb['H13'][i].dict();h22=face_pair(full['records'][i],(2,2));h50=face_pair(full['records'][i],(5,0))
 fs=pull(h13,sigimgs)/beta+GH(sigimgs,siges)+t*pull(h22,sigimgs,True)+t*beta*pull(h50,sigimgs,True)
 fq=pull(h13,quadimgs)/beta+GH(quadimgs,quades)
 add_source('H'+str(i),fs,fq);print('GENERIC BETA H',i,'seconds',time.monotonic()-st,flush=True)
for it,jt in zip(ev['items'],em['items']):
 si=int(it['square_component']);qi=it['Q_index'];ck('full even module source order '+str(qi),qi==jt['Q_index'])
 def even0(imgs,es):return cv(it['square_unit'])*es[si]**2+sum((pull(ff.dict(),imgs)*es[j] for j,ff in enumerate(it['A_p_quadratic'])),T(0))+pull(it['pure_p4'].dict(),imgs)
 # Retain the generically nonzero sigma*RR-linear*e term from the COMPLETE identity.
 extra=sum((pull(ff.dict(),sigimgs)*siges[j] for j,ff in enumerate(jt['B_RR_linear'])),T(0));tail=part(it['full_source'],(2,1));fs=even0(sigimgs,siges)+t*(extra+pull(tail,sigimgs));fq=even0(quadimgs,quades);add_source('Q'+str(qi),fs,fq)
rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'all22_true_sigma_columns_beta_independent':True,'all22_true_R_compensated_kernel_quadratics_beta_independent':True,'original_all_layer_even_extra_retained':True,'all_four_roots':True,'scope':__doc__}
save({'record':rec,'checks':checks,'field':K,'root_ring':E,'labels':labels,'true_R_range':A,'sigma_R_log':siglog,'second_R_log':quadlog,'sigma_column':vector(E,sigma),'kernel_quadratic':vector(E,quadratic),'beta':beta},str(p/'independent_critical_inactive_fullRR_beta_polynomial_reader.sobj'),compress=False);(p/'independent_critical_inactive_fullRR_beta_polynomial_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,flush=True)
