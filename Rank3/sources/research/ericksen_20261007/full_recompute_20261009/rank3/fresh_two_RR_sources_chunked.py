"""Rebuild all 25x5 critical/subcritical two-RR normal matrices, all charts.

Only fresh Section2-derived R/Q/H data supply source coefficients. The three
coordinate substitutions preserve t,w,s and the free even variable e. Every
entry, including kappa and eta columns, is compared only after reconstruction.
"""
from sage.all import *
from pathlib import Path
import time,json,sys,hashlib
out=Path(__file__).absolute().parent;research=out.parent.parent;st=time.monotonic();checks=[]
selected_case,start_row,end_row=map(int,sys.argv[1:4]); assert 0<=selected_case<3 and 0<=start_row<end_row<=18
chunk=out/'two_RR_chunks';chunk.mkdir(exist_ok=True)
inputs=[out/'fresh_global_polynomials.sobj',out/'fresh_global_source_DAG.sobj']+[out/'critical_rows'/('row'+str(j)+'.sobj') for j in range(18)]
fingerprints={str(v.relative_to(out)):hashlib.sha256(v.read_bytes()).hexdigest() for v in inputs}
f=load(str(out/'fresh_global_polynomials.sobj'));dag=load(str(out/'fresh_global_source_DAG.sobj'))
original=f['ring'];old=original.base_ring();K=old;a,b,h,z=K.gens();r=(a-b)/(a*(b-1));q=1+r;H=r*h;Z=r*z/q
P=PolynomialRing(K,names=('t','w','s','e'));t,w,s,e=P.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
def cv(c):return K(c)
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def face(poly,i,j):return original({tuple(ex):c for ex,c in poly.dict().items() if (sum(ex[:6]),sum(ex[6:]))==(i,j)})
pool=[original(v) for v in dag['polynomial_pool']]
Hrows=[]
for j,item in enumerate(dag['records']):
 row=load(str(out/'critical_rows'/('row'+str(j)+'.sobj')))
 Hrows.append(dict(A=row['cofactors'],H41=row['source_faces'][4,1],H13=row['source_faces'][1,3],H50_DAG=item['full_signed_product_DAG']))

# Exact Q12=L cofactor identity, without any saved compatibility object.
Q=f['Q0'];Q12=face(Q,1,2);x=original.gens();cof=[]
for pj,rj in [(0,0),(1,0),(2,1)]:
 dd={}
 for ex,c in Q12.dict().items():
  ee=list(ex)
  if ee[pj]:ee[pj]-=1;ee[6+rj]-=1;ck('Q0 cofactor legal exponents',min(ee)>=0);dd[tuple(ee)]=c
 cof.append(original(dd))
ell=f['L']*vector(original,x[6:]);ck('complete Q12 exact module',sum(cof[j]*ell[j] for j in range(3))==Q12)
DC=f['D']*f['C'];Dd=f['D']*f['offset'];Q21=face(Q,2,1);Q30=face(Q,3,0)
charts=[([w,s,P.one(),P.zero(),t,P.zero(),P.zero(),P(beta)*t,P.one()],2,7,((5,t*w/alpha),(3,s/alpha),(6,t))),([P.one(),P.zero(),w,s,P.zero(),t,P(alpha)*t,P.zero(),P.one()],0,6,((1,s/beta),(4,w*t/beta),(7,t))),([P.zero(),P.one(),P.zero(),t,w,s,P(alpha/beta)*t,P.one(),P.zero()],1,6,((0,beta*s),(2,beta*w*t),(8,t)))]
matrices=[]
for case,(vals,free,rrder,etap) in enumerate(charts):
 if case!=selected_case:continue
 active=[j for j in range(3) if j!=free]
 def pull(poly,derivative=None):
  if derivative is not None:poly=poly.derivative(x[derivative])
  return sum((P(cv(c))*prod(vals[k]**int(ex[k]) for k in range(9)) for ex,c in poly.dict().items() if all(not ex[k] or vals[k]!=0 for k in range(9))),P.zero())
 ck('literal three L vanish case'+str(case),all(pull(z)==0 for z in ell))
 rows=[]
 for j,R in enumerate(f['sources']):
  ck('whole R constant zero'+str((case,j)),pull(DC[j,free])*e+pull(Dd[j])==0)
  eta=sum((pull(DC[j,free],k)*e+pull(Dd[j],k))*coeff for k,coeff in etap)
  rows.append([pull(DC[j,active[0]]),pull(DC[j,active[1]]),P.zero(),pull(face(R,1,2)),eta])
 ck('whole Q0 constant zero'+str(case),pull(cof[free])*e+pull(Q21)==0)
 kap=e*(pull(cof[free],rrder)*e+pull(Q21,rrder))+pull(Q30)
 eta=sum((pull(cof[free],k)*e+pull(Q21,k))*coeff for k,coeff in etap)
 rows.append([pull(cof[active[0]]),pull(cof[active[1]]),kap,P.zero(),eta])
 for j,it in enumerate(Hrows):
  if not start_row<=j<end_row:continue
  checkpoint=chunk/('fresh_case%d_H%d.sobj'%(case,j))
  if checkpoint.is_file():
   saved=load(str(checkpoint));ck('fresh-only resume hashes '+str((case,j)),saved['fingerprints']==fingerprints);rows.append(saved['row']);continue
  print('fresh twoRR original field',case,j,round(time.monotonic()-st,2),flush=True)
  AH,hh=it['A'],it['H41'];ck('whole H constant zero'+str((case,j)),pull(AH[free])*e+pull(hh)==0)
  h50pull=sum((sign*pull(face(pool[li],pd,0))*pull(face(pool[ri],5-pd,0)) for sign,li,ri in it['H50_DAG'] for pd in range(6)),P.zero())
  kap=e*(pull(AH[free],rrder)*e+pull(hh,rrder))+h50pull
  eta=sum((pull(AH[free],k)*e+pull(hh,k))*coeff for k,coeff in etap)
  rows.append([pull(AH[active[0]]),pull(AH[active[1]]),kap,pull(it['H13']),eta])
  save(dict(row=rows[-1],fingerprints=fingerprints),str(checkpoint),compress=False)
 M=matrix(P,rows);matrices.append(M)
 name='independent_subcritical_two_RR_full_normal'+('_case'+str(case) if case<2 else '')+'.sobj'
 archived=load(str(research/'round10'/name));ak=archived['ring'].base_ring();amp=ak.ring().hom([a,r,H,Z],K)
 def ac(c):c=ak(c);return amp(c.numerator())/amp(c.denominator())
 row_ids=list(range(7))+list(range(7+start_row,7+end_row))
 for ii,i in enumerate(row_ids):
  for j in range(5):ck('ALL twoRR source entries'+str((case,i,j)),M[ii,j]==P({tuple(map(int,ex)):ac(c) for ex,c in archived['matrix'][i,j].dict().items()}))
 print('fresh complete twoRR chart ORIGINAL MATERIAL',case,round(time.monotonic()-st,2),flush=True)
name='case%d_rows%d_%d'%(selected_case,start_row,end_row)
save(dict(field=K,ring=P,matrix=matrices[0],row_ids=row_ids,case=selected_case,charts=charts),str(chunk/(name+'.sobj')),compress=False)
inputs=[out/'fresh_global_polynomials.sobj',out/'fresh_global_source_DAG.sobj']+[out/'critical_rows'/('row'+str(j)+'.sobj') for j in range(18)]
(chunk/(name+'.json')).write_text(json.dumps(dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,case=selected_case,rows=row_ids,source_sha256={str(v.relative_to(out)):hashlib.sha256(v.read_bytes()).hexdigest() for v in inputs},source='Fresh Section2-derived polynomial DAG. Exact homogeneous-product convolution H50; original coefficient field preserved; final comparison with archived matrix only after full source differentiation.'),indent=2))
print('PASS',name,len(checks),time.monotonic()-st,flush=True)
