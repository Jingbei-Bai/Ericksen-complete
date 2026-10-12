from sage.all import *
from pathlib import Path
import time,json,os,sys
case=int(sys.argv[1])
p=Path(__file__).resolve().parent;st=time.monotonic();src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));K=src['ring'].base_ring();a,r,H,Z=K.gens();q=1+r;alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));O=src['ring'];x=O.gens();P=PolynomialRing(K,names=('nr_t','nr_w','nr_s','nr_e'));t,w,s,e=P.gens();charts=[([w,s,P.one(),P.zero(),t,P.zero(),P.zero(),P(beta)*t,P.one()],2,7,((5,t*w/alpha),(3,s/alpha),(6,t))),([P.one(),P.zero(),w,s,P.zero(),t,P(alpha)*t,P.zero(),P.one()],0,6,((1,s/beta),(4,w*t/beta),(7,t))),([P.zero(),P.one(),P.zero(),t,w,s,P(alpha/beta)*t,P.one(),P.zero()],1,6,((0,beta*s),(2,beta*w*t),(8,t)))]; vals,free,rrder,etap=charts[case];active=[j for j in range(3) if j!=free];checks=[];labels=[];rows=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
def cp(f):return O({tuple(map(int,ex)):K(c) for ex,c in f.dict().items()})
def face(f,i,j):return O({tuple(map(int,ex)):K(c) for ex,c in f.dict().items() if sum(ex[:6])==i and sum(ex[6:])==j})
def pull(f):
 out=P.zero()
 for ex,c in f.dict().items():
  mon=prod(vals[k]**int(ex[k]) for k in range(9))
  if mon:out+=K(c)*mon
 return out
comp=load(str(p/'independent_delta_Q0_kernel_compatibility.sobj'));off=vector(O,[cp(f) for f in comp['offset']]);C=matrix(O,src['C']);D=matrix(O,src['D']);DC=D*C;Dd=D*off
ck('exact two-RR L kernel',all(pull(f)==0 for f in matrix(O,src['L'])*vector(O,x[6:])))
for j,Rf in enumerate(src['full_R_sources']):
 Rf=cp(Rf);zero=pull(DC[j,free])*e+pull(Dd[j]);ck('R boundary zero '+str(j),zero==0)
 eta=sum((pull(DC[j,free].derivative(x[k]))*e+pull(Dd[j].derivative(x[k])))*coeff for k,coeff in etap)
 rows.append([pull(DC[j,active[0]]),pull(DC[j,active[1]]),P.zero(),pull(face(Rf,1,2)),eta]);labels.append(('R',j))
Q=cp(src['full_Q0']);B=list(map(cp,src['Q0_L_cofactors']));Q21=face(Q,2,1);Q30=face(Q,3,0)
ck('Q0 boundary zero',pull(B[free])*e+pull(Q21)==0)
kappa=e*(pull(B[free].derivative(x[rrder]))*e+pull(Q21.derivative(x[rrder])))+pull(Q30)
eta=sum((pull(B[free].derivative(x[k]))*e+pull(Q21.derivative(x[k])))*coeff for k,coeff in etap)
rows.append([pull(B[active[0]]),pull(B[active[1]]),kappa,P.zero(),eta]);labels.append(('Q0',0))
mod=load(str(p/'independent_delta_H_global_L_module.sobj'));old=mod['original_field'];Hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));pool=full['polynomial_pool'];ph=old.ring().hom([a,a*q/(1+a*r),H/r,q*Z/r],K);cache={}
def cv(c):
 c=old(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
def oldface(dag,i,j):
 dd={}
 for sign,lref,rref in dag:
  for ex,c in pool[lref].items():
   for ey,b in pool[rref].items():
    ee=tuple(int(v)+int(z) for v,z in zip(ex,ey))
    if sum(ee[:6])==i and sum(ee[6:])==j:dd[ee]=dd.get(ee,old.zero())+sign*c*b
 return {ex:c for ex,c in dd.items() if c}
def opull(dd,derivative=None):
 out=P.zero()
 for ex,c in dd.items():
  ex=list(map(int,ex));mul=1
  if derivative is not None:
   mul=ex[derivative]
   if not mul:continue
   ex[derivative]-=1
  mon=prod(vals[k]**ex[k] for k in range(9))
  if mon:out+=mul*cv(c)*mon
 return out
previous_path=p/('independent_subcritical_two_RR_full_normal_case'+str(case)+'.sobj')
if previous_path.exists():
 prior=load(str(previous_path));assert prior['labels'][:7]==labels
 assert list(prior['matrix'][:7,:].rows())==list(matrix(P,rows).rows())
 rows=[list(row) for row in prior['matrix'].rows()];labels=list(prior['labels']);checks=list(prior['checks'])
for j,it in enumerate(Hs['items']):
 if j<len(rows)-7:continue
 AH=it['L_cofactors_original'];hh=it['H41_original'];h50=oldface(full['records'][j]['full_signed_product_DAG'],5,0);h13=oldface(full['records'][j]['full_signed_product_DAG'],1,3)
 ck('H complete boundary zero '+str(j),opull(AH[free])*e+opull(hh)==0)
 kap=e*(opull(AH[free],rrder)*e+opull(hh,rrder))+opull(h50)
 eta=sum((opull(AH[free],k)*e+opull(hh,k))*coeff for k,coeff in etap)
 rows.append([opull(AH[active[0]]),opull(AH[active[1]]),kap,opull(h13),eta]);labels.append(('H',j))
 print('ROW',j,'entryterms',[len(z.dict()) for z in rows[-1]],'seconds',time.monotonic()-st,flush=True)
 tmp=p/('independent_subcritical_two_RR_full_normal_case'+str(case)+'.tmp.sobj');save({'coefficient_field':K,'ring':P,'variables':P.gens(),'p_images':vals[:6],'RR_images':vals[6:],'matrix':matrix(P,rows),'labels':labels,'checks':checks,'column_order':['e'+str(active[0]),'e'+str(active[1]),'kappa','delta','t_times_eta'],'missing_RR':case,'free_e':free,'active_e':active,'eta_tangent':etap,'rr_kappa_index':rrder,'source_files':['independent_delta_subcritical_Q0_full_source.sobj','independent_delta_subcritical_H_initials.sobj','independent_ratio_global_R_derivative_degreefall.sobj'],'scope':'Complete regular R,Q0,H first normal at kappa=delta=eta=e0=e2=0, on the rankD3 missing-RR3 chart. Free e1 and both w,s retained; t nonzero. No inverse C and no w/s division. Additional source even equations and lower-rank loci retained.'},str(tmp),compress=False);os.replace(tmp,p/('independent_subcritical_two_RR_full_normal_case'+str(case)+'.sobj'))
(p/('independent_subcritical_two_RR_full_normal_case'+str(case)+'.json')).write_text(json.dumps({'checks':len(checks),'rows':len(rows),'seconds':time.monotonic()-st,'all_sources_complete':len(rows)==25},indent=2)+'\n');print('DONE',len(rows),len(checks),time.monotonic()-st,flush=True)
