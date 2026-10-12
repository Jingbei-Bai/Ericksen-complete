from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks={}
def ck(n,v):
 if not v:raise AssertionError(n)
 checks[n]=True
raw=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));d=load(str(p/'independent_ratio_RR_canonical_tangent_obstruction.sobj'));K=raw['coefficient_field'];a,r=raw['parameters'];q=1+r;guard=raw['guard'];M=raw['matrix'];rr=M.base_ring().gens();T=PolynomialRing(K,names=('auditX','auditY'));X,Y=T.gens();scales=[rr[2],rr[1],rr[2],rr[0],rr[1],rr[0],rr[0]*rr[1]*rr[2]];rows=[]
for i,row in enumerate(M.rows()):
 fs=[f*s for f,s in zip(row,scales)];mons=[tuple(map(int,e)) for f in fs for e in f.dict()];pa=tuple(v%2 for v in mons[0]);ck(f'parity_{i}',all(tuple(v%2 for v in e)==pa for e in mons))
 out=[]
 for f in fs:
  poly=T.zero()
  for e,c in f.dict().items():
   ex=tuple((int(ej)-pj)//2 for ej,pj in zip(e,pa));ck(f'exponent_nonnegative_{i}_{len(checks)}',min(ex)>=0);poly+=c*(q/r)**ex[0]*q**ex[1]*X**ex[0]*Y**ex[1]
  out.append(poly)
 rows.append(out)
N=matrix(T,rows);saved=d['torus_matrix'];ck('torus_matrix_source',all(N[i,j]==T({tuple(map(int,e)):K(c) for e,c in saved[i,j].dict().items()}) for i in range(28) for j in range(7)))
N0=matrix(K,28,7,[f(-1,-1) for f in N.list()]);ck('constant_source',N0==d['constant_matrix']);kr=d['kernel_rows'];ck('whole_kernel',N0*kr.transpose()==0 and kr.rank()==3)
cols=list(d['pivot_columns']);ids=list(d['pivot_rows']);C=N0.matrix_from_columns(cols);P=C.matrix_from_rows(ids);det=P.det();ck('pivot_source',det==d['pivot_determinant']);ck('pivot_protected',all(guard.quo_rem(f)[1]==0 for s in [det.numerator(),det.denominator()] for f,e in s.factor()))
Bi=P.inverse();B=matrix(K,28,6,0)
for direction,var in enumerate((X,Y)):
 Nd=matrix(K,28,7,[f.derivative(var)(-1,-1) for f in N.list()])
 for j,v in enumerate(kr.rows()):
  rawv=Nd*v;correction=Bi*vector(K,[rawv[i] for i in ids]);remainder=rawv-C*correction;B.set_column(3*direction+j,remainder)
  ck(f'projected_pivot_zero_{direction}_{j}',all(remainder[i]==0 for i in ids))
ck('obstruction_source',B==d['obstruction_matrix']);rs=list(d['rank_minor_rows']);dd=B.matrix_from_rows(rs).det();ck('six_minor_source',dd==d['rank_minor'] and dd!=0);ck('six_minor_protected',all(guard.quo_rem(f)[1]==0 for s in [dd.numerator(),dd.denominator()] for f,e in s.factor()));ck('short_minor_identity',dd==QQ(2048)/729*q**6/r**5)
# Verify the old seven complete sources meet the exact normal-order and
# parity conditions used in the formal-arc inference.
f=load(str(p/'independent_ratio_weighted_degreefall.sobj'));ss=[f['Q0']]+[i['full'] for i in f['items']];wt=(1,)*6+(2,)*3
for j,g in enumerate(ss):
 ck(f'full_weight5_{j}',all(sum(int(ei)*wi for ei,wi in zip(e,wt))<=5 for e in g.dict()))
 ck(f'RRdegree2_{j}',all(sum(map(int,e[6:]))<=2 for e in g.dict()))
 ck(f'normal_odd_{j}',all((5-2*sum(map(int,e[6:])))%2==1 and 5-2*sum(map(int,e[6:]))>=1 for e in g.dict()))
report={'status':'PASS','count':len(checks),'checks':checks,'six_minor':str(dd.factor()),'seconds':time.monotonic()-st,'scope':'Original28 source torus matrix, full kernel, protected4pivot, exact derivative/projection and protected6minor replayed. At FIXED cyclic base the complete odd normal sources imply RR motion valuation at least twice leading kernel valuation. No scalar removal of all kernel corrections, no affine-lift exclusion.'}
(p/'independent_RR_tangent_obstruction_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks,'six_minor':dd},str(p/'independent_RR_tangent_obstruction_reader.sobj'));print({k:v for k,v in report.items() if k!='checks'},flush=True)
