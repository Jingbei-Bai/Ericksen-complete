"""Fresh complete primary normal system from actual joint fields and R sources."""
from sage.all import *
from pathlib import Path
import json,time,hashlib
out=Path(__file__).absolute().parent;st=time.monotonic();checks=[]
f=load(str(out/'fresh_global_polynomials.sobj'));d=load(str(out/'fresh_joint_fields.sobj'))
R=f['ring'];K=R.base_ring();a,b,h,z=K.gens();xx=R.gens();sources=[f['Q0']]+list(f['sources']);hat=f['fields'];base=f['base_fields'];wt=(1,)*6+(2,)*3
T=PolynomialRing(K,names=('auditRR1','auditRR2','auditRR3'));rr=T.gens()
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def normal(f,degree):
 row=[T.zero() for j in range(7)]
 for e,c in f.dict().items():
  if sum(e[6:])!=degree:continue
  ck('normal p degree '+str(len(checks)),sum(e[:6])<=1)
  k=next((k for k in range(6) if e[k]),6);row[k]+=T({tuple(e[6:]):c})
 return row
def pure(f,degree):return T({tuple(e[6:]):c for e,c in f.dict().items() if not any(e[:6]) and sum(e[6:])==degree})
N0=matrix(T,[normal(f,2) for f in sources]);L=matrix(T,[[pure(v,1) for v in row] for row in base.rows()]);RR=matrix(T,[[pure(hat[i,6+j],2) for j in range(3)] for i in range(3)])
DP=[matrix(T,[normal(d['fields'][i,j],1) for j in range(6)]) for i in range(3)]
for i in range(3):
 for j in range(9):ck('full gauged weight bound '+str((i,j)),all(sum(e[k]*wt[k] for k in range(9))<=(3 if j<6 else 4) for e in hat[i,j].dict()))
def bd(g,v):return T({tuple(e):c.derivative(v) for e,c in g.dict().items() if c.derivative(v)})
def deriv(i,g):return L[i,0]*bd(g,h)+L[i,1]*bd(g,z)+sum(RR[i,j]*g.derivative(rr[j]) for j in range(3))
rows=list(N0.rows())
for j in range(7):
 for i in range(3):
  row=[deriv(i,N0[j,k])+sum(N0[j,l]*DP[i][l,k] for l in range(6)) for k in range(7)];rows.append(vector(T,row))
  print('fresh higher normal',j,i,round(time.monotonic()-st,2),flush=True)
M=matrix(T,rows)
ck('seven original normal rows',N0.dimensions()==(7,7))
ck('all twenty-one ordered first derivatives',M.dimensions()==(28,7))
ck('complete Q0 is nonzero',bool(f['Q0']))
qex=sorted(f['Q0'].dict())[0];qcoeff=f['Q0'].dict()[qex];gauges=[]
for i in range(3):
 difference=R(d['fields'][i,6+i])-hat[i,6+i]
 gauge=K(difference.dict().get(qex,K.zero())/qcoeff)
 ck('whole original Q0 gauge '+str(i),difference==gauge*f['Q0']);gauges.append(gauge)
save(dict(ring=T,parameters=(a,b,h,z),old_normal=N0,normal_matrix=M,base_fields=L,RR_fields=RR,p_derivative_blocks=DP),str(out/'fresh_primary_normals.sobj'),compress=False)
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,source='Fresh actual joint fields and complete R/Q0 sources, seven original normal rows and all twenty-one ordered first derivatives.',all28rows=True,Q0_gauge_coefficients=[str(v) for v in gauges],inputs=[dict(path=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in [out/'fresh_global_polynomials.sobj',out/'fresh_joint_fields.sobj']])
(out/'fresh_primary_normal_sources.json').write_text(json.dumps(report,indent=2));print('PASS',len(checks),report['seconds'],flush=True)
