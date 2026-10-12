from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks={};reports=[]
def ck(n,b):
 if not b:raise AssertionError(n)
 checks[n]=True
raw=load(str(p/'joint_fields_direct.sobj'));fd=load(str(p/'independent_ratio_weighted_degreefall.sobj'));kd=load(str(p/'independent_ratio_diagonal_Q0_reduction.sobj'));sec=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));oldK=raw['ring'].base_ring();a,b,h,z=raw['parameters'];K=sec['coefficient_field'];aa,r=sec['parameters'];q=1+r;hp=oldK.ring().hom(sec['old_parameter_map'],K);w=(1,)*6+(2,)*3
J=PolynomialRing(K,names=('audit_dh','audit_dz'));dh,dz=J.gens();zero=(0,)*9;valcache={};jet1cache={};jet2cache={}
def ev(c):
 c=oldK(c)
 if c not in valcache:valcache[c]=hp(c.numerator())/hp(c.denominator())
 return valcache[c]
def jt1(c):
 c=oldK(c)
 if c not in jet1cache:jet1cache[c]=J(ev(c))+ev(c.derivative(h))*dh+ev(c.derivative(z))*dz
 return jet1cache[c]
def jt2(c):
 c=oldK(c)
 if c not in jet2cache:
  ch=c.derivative(h);cz=c.derivative(z);jet2cache[c]=jt1(c)+ev(ch.derivative(h))/2*dh**2+ev(ch.derivative(z))*dh*dz+ev(cz.derivative(z))/2*dz**2
 return jet2cache[c]
def trunc1(f):return J({tuple(map(int,e)):c for e,c in f.dict().items() if sum(e)<=1 and c})
def wd(e):return sum(w[j]*e[j] for j in range(9))
def tod(f):return {tuple(map(int,e)):oldK(c) for e,c in f.dict().items() if c}
def put(dic,e,vv):
 if vv:dic[e]=dic.get(e,J.zero())+vv
Q=tod(fd['Q0']);fields=[];bases=[]
for i in range(3):
 row=[]
 for j in range(9):
  vals=tod(raw['fields'][i,j])
  if j==6+i:
   for e,c in Q.items():vals[e]=vals.get(e,oldK.zero())-kd['items'][i]['kappa']*c
  vals={e:c for e,c in vals.items() if c};ck('entire_field_weight_bound_%d_%d'%(i,j),all(wd(e)-w[j]<=2 for e in vals))
  for e in vals:ck('native_tuple_field_%d_%d_%s'%(i,j,e),type(e)==tuple and all(type(k)==int for k in e))
  row.append({e:jt1(c) for e,c in vals.items() if wd(e)-w[j] in (1,2)})
 fields.append(row);bb=[]
 for j in range(2):
  vals=tod(raw['base_fields'][i,j]);ck('entire_base_weight_bound_%d_%d'%(i,j),all(wd(e)<=2 for e in vals));bb.append({e:jt1(c) for e,c in vals.items() if wd(e) in (1,2)})
 bases.append(bb)
print('complete field filtration and independent Taylor coefficient jets prepared',round(time.monotonic()-st,2),flush=True)
# Two actions are performed in a truncated Taylor polynomial algebra in dh,dz.
# This independently produces the differentiated coefficients rather than copying
# the creator's explicit three-/six-tuple product-rule implementation.
def act1(source,i):
 out={}
 for e,c in source.items():
  c2=jt2(c)
  for k in range(9):
   if not e[k]:continue
   ee=list(e);ee[k]-=1
   for f,f1 in fields[i][k].items():
    nu=tuple(ee[j]+f[j] for j in range(9))
    if wd(nu) in (6,7):put(out,nu,e[k]*trunc1(c2*f1))
  for k,dvar in enumerate((dh,dz)):
   der=c2.derivative(dvar)
   for f,f1 in bases[i][k].items():
    nu=tuple(e[j]+f[j] for j in range(9))
    if wd(nu) in (6,7):put(out,nu,trunc1(der*f1))
 return {e:trunc1(c) for e,c in out.items() if c}
def act2(source,i):
 out={}
 def add(e,c):out[e]=out.get(e,K.zero())+c
 for e,c1 in source.items():
  c0=c1.constant_coefficient()
  for k in range(9):
   if not e[k]:continue
   ee=list(e);ee[k]-=1
   for f,f1 in fields[i][k].items():
    nu=tuple(ee[j]+f[j] for j in range(9))
    if wd(nu) in (8,9):add(nu,e[k]*c0*f1.constant_coefficient())
  for k,dvar in enumerate((dh,dz)):
   dc=c1.monomial_coefficient(dvar)
   for f,f1 in bases[i][k].items():
    nu=tuple(e[j]+f[j] for j in range(9))
    if wd(nu) in (8,9):add(nu,dc*f1.constant_coefficient())
 return {e:c for e,c in out.items() if c}
A=PolynomialRing(K,names=('audit_t0','audit_t1','audit_t2'));tt=A.gens();sq1=-q/r;sq2=-q
# Independent algebraic-root reduction. Introduce roots s1,s2,c with
# s1^2=-q/r, s2^2=-q, c^2=s1*s2, RR3=1. Do not set c or J=s1*s2 to one.
def radical_coefficient(co,es1,es2,ec,n):
 k,rc=divmod(int(ec),2);es1+=k;es2+=k;k1,r1=divmod(int(es1),2);k2,r2=divmod(int(es2),2)
 ck(n+'_root_basis_rational',(r1,r2,rc)==(0,0,0));return co*sq1**k1*sq2**k2
normal_count=0
def independent_normal(top,low,parity,label):
 out=[A.zero() for _ in range(6)];cone=A.zero()
 for e,co in top.items():
  n=sum(e[:6]);ck(label+'_odd_p_top_'+str(e),n%2==1);rho=(e[6]+e[3]+e[5],e[7]+e[1]+e[4],e[8]+e[0]+e[2]);ke=(e[0]+e[5],e[1]+e[3],e[2]+e[4]);d0=rho[0]-parity[0];d1=rho[1]-parity[1]
  val=radical_coefficient(co,d0,d1,1-n,label+'_top_'+str(e));mon=A.monomial(*ke);cone+=val*mon
  for j in range(3):
   if e[j]:kk=list(ke);kk[j]-=1;out[j]+=e[j]*val*A.monomial(*kk)
  # RR1=s1*(1-Xnormal/2)+O(Xnormal^2), and analogous RR2.
  # The parity row unit is differentiated as well; c is fixed at the cone.
  out[4]+=-K(d0)/2*val*mon;out[5]+=-K(d1)/2*val*mon
 for e,co in low.items():
  n=sum(e[:6]);ck(label+'_even_p_tau_'+str(e),n%2==0);rho=(e[6]+e[3]+e[5],e[7]+e[1]+e[4],e[8]+e[0]+e[2]);ke=(e[0]+e[5],e[1]+e[3],e[2]+e[4]);val=radical_coefficient(co,rho[0]-parity[0],rho[1]-parity[1],2-n,label+'_tau_'+str(e));out[3]+=val*A.monomial(*ke)
 ck(label+'_whole_cone_zero',cone==0);return out
sources=[fd['Q0']]+[it['full'] for it in fd['items']]
allrows=[]
for sn,f in enumerate(sources):
 saved=load(str(p/('independent_ratio_canonical_second_prolongation_source%d.sobj'%sn)));ck('source_complete_nine_'+str(sn),len(saved['items'])==9 and len(saved['first_source_jets'])==3);ss=tod(f);ck('entire_source_weight_bound_'+str(sn),all(wd(e)<=5 for e in ss));selected={e:c for e,c in ss.items() if wd(e) in (4,5)};first=[]
 for inner in range(3):
  jj=act1(selected,inner);first.append(jj);old={tuple(map(int,e)):tuple(K(v) for v in vals) for e,vals in saved['first_source_jets'][inner].items()};expected={e:(cc.constant_coefficient(),cc.monomial_coefficient(dh),cc.monomial_coefficient(dz)) for e,cc in jj.items()};ck('first_entire_exponent_support_%d_%d'%(sn,inner),set(old)==set(expected))
  for e,vals in expected.items():
   for k in range(3):ck('first_independent_Taylor_coeff_%d_%d_%s_%d'%(sn,inner,e,k),vals[k]==old[e][k])
  for outer in range(3):
   dd=act2(jj,outer);top={e:c for e,c in dd.items() if wd(e)==9};low={e:c for e,c in dd.items() if wd(e)==8};it=saved['items'][3*inner+outer];label='source%d_inner%d_outer%d'%(sn,inner,outer);ck(label+'_indices',it['source_index']==sn and it['inner']==inner and it['outer']==outer)
   for name,new in [('top',top),('tau',low)]:
    old={tuple(map(int,e)):K(c) for e,c in it[name].items() if c};ck(label+'_'+name+'_complete_support',set(new)==set(old))
    for e,c in new.items():ck(label+'_'+name+'_original_coeff_'+str(e),c==old[e])
   nr=independent_normal(top,low,tuple(map(int,it['parity'])),label)
   for j in range(6):
    old=A({tuple(map(int,e)):K(c) for e,c in it['normal_row'][j].dict().items()});ck(label+'_normal_column_'+str(j),nr[j]==old)
   allrows.append(nr);reports.append({'source':sn,'inner':inner,'outer':outer,'top_terms':len(top),'tau_terms':len(low),'normal_replayed':True})
  print('source',sn,'inner',inner,'all three genuine second actions and radical-normal columns PASS',round(time.monotonic()-st,2),flush=True)
 interim={'status':'PARTIAL_EXACT_REPLAY','completed_sources':sn+1,'completed_rows':len(allrows),'count':len(checks),'seconds':time.monotonic()-st};(p/'independent_second_prolongation_reader_progress.json').write_text(json.dumps(interim,indent=2)+'\n')
ck('complete_63_rows',len(allrows)==63)
report={'status':'PASS','count':len(checks),'seconds':time.monotonic()-st,'records':reports,'scope':'All63 genuine hatFi hatFj sources independently reconstructed via truncated Taylor polynomial algebra before cyclic substitution, with full weighted-filtration completeness verified from untruncated fields and sources. Weight9 tops/weight8 tau coefficients compared individually; six-column normal conversion independently reduced in roots s1^2=-q/r,s2^2=-q,c^2=s1*s2 with c kept nontrivial. All lower layers capable of contributing and two h,z coefficient derivatives included. Cyclic BASE only, all protected material spectra; necessary normal sources, no rank/fiber exclusion asserted.'}
(p/'independent_second_prolongation_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks,'records':reports,'amplitude_ring':A,'normal_matrix':matrix(A,allrows)},str(p/'independent_second_prolongation_reader.sobj'));print(report['status'],report['count'],report['seconds'],flush=True)
