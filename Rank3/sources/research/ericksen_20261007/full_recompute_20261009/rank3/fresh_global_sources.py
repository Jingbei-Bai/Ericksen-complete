"""Regenerate R0..R5, Q0 and all 18 complete necessary derivative DAGs.

Starts at freshly rebuilt actual joint-coordinate fields and all 21 quartics.
All canceled layers, coefficients, partials and degree bounds are recomputed.
The six source multipliers and three field gauges are solved by coefficient
comparison, not imported. Final historical DAG equality is a cross-check only.
"""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;research=out.parent.parent;start=time.monotonic();checks=[]
d=load(str(out/'fresh_joint_fields.sobj'));R=d['ring'];K=R.base_ring();a,b,h,z=K.gens();xx=R.gens();qs=d['quartic_source_rows'];Q0=qs[0]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def layer(f,degree):return R({tuple(e):c for e,c in f.dict().items() if sum(e)==degree})
def face(f,pd,rd):return R({tuple(e):c for e,c in f.dict().items() if sum(e[:6])==pd and sum(e[6:])==rd})
def weight(f):return max(sum(e[:6])+2*sum(e[6:]) for e in f.dict())
def topweight(f):wd=weight(f);return R({tuple(e):c for e,c in f.dict().items() if sum(e[:6])+2*sum(e[6:])==wd})
def deg(f):return max(map(sum,f.dict()),default=-1)
def bd(f,param):return R({tuple(e):c.derivative(param) for e,c in f.dict().items() if c.derivative(param)})
def prog(n,**kw):print(n,round(time.monotonic()-start,2),kw,flush=True)
sources=[];source_mult=[]
for index in [10,11,13,15,18,19]:
 lead=topweight(qs[index]);qtop=topweight(Q0);found=[]
 for j in range(6):
  test=xx[j]*qtop
  if set(test.dict())!=set(lead.dict()):continue
  ex=next(iter(test.dict()));ratio=lead.dict()[ex]/test.dict()[ex]
  if lead==ratio*test:found.append((j,ratio))
 ck('unique exact degree-fall multiplier '+str(index),len(found)==1)
 j,ratio=found[0];full=qs[index]-ratio*xx[j]*Q0
 ck('degree-fall identity '+str(index),full+ratio*xx[j]*Q0==qs[index]);ck('full R ordinary degree4 '+str(index),deg(full)<=4)
 sources.append(full);source_mult.append((j,ratio));prog('R regenerated',index=index)
r=(a-b)/(a*(b-1));q=1+r;H=r*h;Z=r*z/q;alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
L=matrix(R,[[xx[0],0,-alpha*xx[5]],[xx[1],-alpha/beta*xx[3],0],[0,xx[2],-beta*xx[4]]]);DL=L.det();adjL=L.adjugate()
ck('L determinant',DL==alpha*(xx[0]*xx[3]*xx[4]-xx[1]*xx[2]*xx[5]))
BR=matrix(R,[[face(f,3,1).derivative(x) for x in xx[6:]] for f in sources]);B=[]
def exact_divide(num,den):
 work={tuple(e):c for e,c in num.dict().items()};gd={tuple(e):c for e,c in den.dict().items()};lead=max(gd);lc=gd[lead];qq={};rem={}
 while work:
  e=max(work);c=work[e]
  if all(e[i]>=lead[i] for i in range(9)):
   qe=tuple(e[i]-lead[i] for i in range(9));qc=c/lc;qq[qe]=qq.get(qe,K.zero())+qc
   for ee,cc in gd.items():
    key=tuple(ee[i]+qe[i] for i in range(9));val=work.get(key,K.zero())-qc*cc
    if val:work[key]=val
    else:work.pop(key,None)
  else:rem[e]=work.pop(e)
 return R(qq),R(rem)
for row in BR.rows():
 nums=row*adjL;qq=[]
 for num in nums:
  quotient,rem=exact_divide(num,DL);ck('exact module L division',rem==0);qq.append(quotient)
 B.append(qq)
B=matrix(R,B)
D=matrix(R,6,3,0);D[0,0]=xx[4];D[1,0]=-xx[5];D[2,1]=xx[2];D[5,1]=-xx[3];D[3,2]=xx[0];D[4,2]=xx[1]
def divide_var(f,j):
 ans={}
 for ee,c in f.dict().items():
  ee=list(ee);ck('exact coordinate division',ee[j]>0);ee[j]-=1;ans[tuple(ee)]=c
 return R(ans)
C=matrix(R,[[divide_var(f,j) for f in B.row(ii)] for ii,j in [(0,4),(2,2),(3,0)]])
offset=vector(R,[divide_var(face(sources[ii],4,0),j) for ii,j in [(0,4),(2,2),(3,0)]])
ck('full row factor B=D C',D*C==B);ck('full p4 offset factor',D*offset==vector(R,[face(f,4,0) for f in sources]))
Zsrc=C*L*vector(R,xx[6:])+offset
for j in range(6):ck('Rtop=D Z '+str(j),layer(sources[j],4)==(D*Zsrc)[j])
hat=matrix(R,d['fields']);base=matrix(R,d['base_fields']);gauges=[]
for i in range(3):
 qt=face(Q0,1,2);ft=face(hat[i,6+i],1,2);ex=next(iter(qt.dict()));u=ft.dict()[ex]/qt.dict()[ex]
 ck('gauge whole 1,2 face '+str(i),ft==u*qt);hat[i,6+i]-=u*Q0;gauges.append(u)
units=(2*(H-Z)/(H*Z),-2*(Z-1)/Z,r*(H+1)/(H*q));sig=(2,1,0)
for i in range(3):
 ck('unique cubic diagonal '+str(i),layer(hat[i,6+i],3)==units[i]*Zsrc[sig[i]])
 for j in range(9):ck('gauged field degree '+str((i,j)),deg(hat[i,j])<=(3 if j==6+i else 2))
ecoef=[]
for i in range(3):
 rows=[]
 for j in range(3):
  vv=[R.zero() for _ in range(6)]
  if i==0:vv[3]=C[j,0];vv[4]=C[j,1]
  if i==1:vv[5]=alpha/beta*C[j,1];vv[2]=C[j,2]
  if i==2:vv[1]=alpha*C[j,0];vv[0]=-beta*C[j,2]
  ck('exact E leading module '+str((i,j)),sum(vv[k]*D[k,sig[i]] for k in range(6))==(C*L.column(i))[j]);rows.append(vv)
 ecoef.append(rows)
partials=[[f.derivative(v) for v in xx] for f in sources];basepartials=[[bd(f,v) for v in [h,z]] for f in sources]
pool=[];refs={}
def ref(f):
 key=f
 if key not in refs:refs[key]=len(pool);pool.append({tuple(map(int,e)):c for e,c in f.dict().items()})
 return refs[key]
def rr(values):return [ref(f) for f in values]
records=[]
for i in range(3):
 for j in range(6):
  correction=[sum(D[j,k]*ecoef[i][k][ell] for k in range(3)) for ell in range(6)]
  dag=[(1,hat[i,v],partials[j][v]) for v in range(9) if partials[j][v]]+[(1,base[i,v],basepartials[j][v]) for v in range(2) if basepartials[j][v]]+[(-1,units[i]*correction[ell],sources[ell]) for ell in range(6) if correction[ell]]
  ck('complete degree6 cancellation '+str((i,j)),layer(sources[j],4).derivative(xx[6+i])==sum(D[j,k]*(C*L.column(i))[k] for k in range(3)))
  records.append(dict(direction=i,source_position=j,full_signed_product_DAG=[(sgn,ref(left),ref(right)) for sgn,left,right in dag],degree_bound=5))
payload=dict(coefficient_field=K,original_parameters=(a,b,h,z),variable_names=R.variable_names(),polynomial_pool=pool,records=records,full_source_refs=rr(sources),source_Q0=Q0,full_first_gauged_field_refs=[rr(row) for row in hat.rows()],full_base_field_refs=[rr(row) for row in base.rows()],D_refs=[rr(row) for row in D.rows()],C_refs=[rr(row) for row in C.rows()],L_refs=[rr(row) for row in L.rows()],units=units,sigma=sig,source_order=(10,11,13,15,18,19),checks=checks,
 provenance='Regenerated directly from fresh Section 2 state, printed RR formulas and exact coefficient solving. No historical source arrays were used.')
save(payload,str(out/'fresh_global_source_DAG.sobj'),compress=False)
save(dict(ring=R,parameters=(a,b,h,z),Q0=Q0,sources=sources,fields=hat,base_fields=base,D=D,C=C,L=L,offset=offset),str(out/'fresh_global_polynomials.sobj'),compress=False)
record=dict(status='PASS',count=len(checks),seconds=time.monotonic()-start,checks=checks,source_count=6,derivative_DAG_count=18,pool_count=len(pool),no_material_specialization=True,imported_historical_source_arrays=0)
(out/'fresh_global_sources.json').write_text(json.dumps(record,indent=2));prog('PASS',count=len(checks))
