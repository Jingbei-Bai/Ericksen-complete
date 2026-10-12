from sage.all import *
from itertools import combinations
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks={};payload=[]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks[n]=True
raw=load(str(p/'independent_ratio_RR_cyclic_full_matrix.sobj'));M=raw['matrix'];K=raw['coefficient_field'];a,r=raw['parameters'];q=1+r;RR=M.base_ring();rr=RR.gens();guard=raw['guard'];Q=PolynomialRing(QQ,names=('auditW1','auditW2'));x,y=Q.gens()
for kind,filename in [('p','independent_ratio_RR_cyclic_pblock_minors.sobj'),('aug','independent_ratio_RR_cyclic_block_minors.sobj')]:
 d=load(str(p/filename))
 for bi,it in enumerate(d['items']):
  ids=list(it['source_rows']);cols=list(it['source_columns']);records={tuple(i['selection']):i for i in it['records']};seen=[]
  for sel in combinations(range(8),len(cols)):
   mat=M.matrix_from_rows_and_columns([ids[j] for j in sel],cols);f=mat.det()
   if not f:ck(f'{kind}_{bi}_{sel}_zero_not_saved',sel not in records);continue
   ck(f'{kind}_{bi}_{sel}_has_record',sel in records);rec=records[sel];poly=RR.zero()
   for e,c in rec['polynomial'].dict().items():poly+=K(c)*(r/q*rr[0]**2)**int(e[0])*(rr[1]**2/q)**int(e[1])*rr[2]**(2*int(e[2]))
   rebuilt=rec['parameter_scale']*prod(rr[j]**int(rec['RR_monomial'][j]) for j in range(3))*poly;ck(f'{kind}_{bi}_{sel}_full_source',f==rebuilt)
   for side in (rec['parameter_scale'].numerator(),rec['parameter_scale'].denominator()):
    for fac,e in side.factor():ck(f'{kind}_{bi}_{sel}_protected_'+str(fac),guard.quo_rem(fac)[1]==0)
   seen.append(rec['polynomial'])
  ck(f'{kind}_{bi}_all_saved_covered',len(seen)==len(records));ck(f'{kind}_{bi}_polynomial_set',set(seen)==set(it['polynomials']))
  fs=[Q({(int(e[0]),int(e[1])):QQ(c) for e,c in f.dict().items()}) for f in it['polynomials']]
  for j,f in enumerate(fs):ck(f'{kind}_{bi}_{j}_reverse_zero',f(-1,-1)==0)
  # A small rational linear system gives explicit source combinations for the
  # reduced-point generators. No Groebner result is accepted as membership.
  solutions=None
  for degree in range(max(f.total_degree() for f in fs),7):
   entries=[]
   for j,f in enumerate(fs):
    for u in range(degree-f.total_degree()+1):
     for v in range(degree-f.total_degree()-u+1):entries.append((j,u,v,f*x**u*y**v))
   mon=sorted(set(tuple(map(int,e)) for _,_,_,f in entries for e in f.dict())|{(0,0),(1,0),(0,1)})
   A=matrix(QQ,[[f.monomial_coefficient(x**e[0]*y**e[1]) for _,_,_,f in entries] for e in mon]);targets=[vector(QQ,[g.monomial_coefficient(x**e[0]*y**e[1]) for e in mon]) for g in (x+1,y+1)]
   try:cofs=[A.solve_right(t) for t in targets]
   except ValueError:continue
   identities=[]
   for k,g in enumerate((x+1,y+1)):
    cs=[Q.zero() for _ in fs]
    for coeff,(j,u,v,f) in zip(cofs[k],entries):cs[j]+=coeff*x**u*y**v
    ck(f'{kind}_{bi}_forward_identity_{k}',sum((c*f for c,f in zip(cs,fs)),Q.zero())==g);identities.append(cs)
   solutions=identities;break
  ck(f'{kind}_{bi}_forward_found',solutions is not None)
  reverse=[]
  for j,f in enumerate(fs):
   a0,re=(f-Q(f(-1,y))).quo_rem(x+1);ck(f'{kind}_{bi}_reverse_div1_{j}',re==0);b0,re=Q(f(-1,y)).quo_rem(y+1);ck(f'{kind}_{bi}_reverse_div2_{j}',re==0);ck(f'{kind}_{bi}_reverse_identity_{j}',f==a0*(x+1)+b0*(y+1));reverse.append((a0,b0))
  payload.append({'kind':kind,'block':bi,'sources':fs,'forward_multipliers':solutions,'reverse_multipliers':reverse,'degree':degree,'original_source_file':filename});print(kind,bi,'certificate degree',degree,'seconds',round(time.monotonic()-st,2),flush=True)
report={'status':'PASS','count':len(checks),'checks':checks,'seconds':time.monotonic()-st,'scope':'Every original p-block and augmented-block maximal minor reconstructed with only RR-torus and old spectral units; all six torus ideals exactly equal (W1+1,W2+1) by explicit forward and reverse source identities. This supplies scheme-theoretic reducedness, including the union of all augmented minors. Coordinate-plane coverage is a separate previously audited certificate.'}
(p/'independent_RR_special_fiber_scheme_reader.json').write_text(json.dumps(report,indent=2)+'\n');save({'checks':checks,'items':payload},str(p/'independent_RR_special_fiber_scheme_reader.sobj'));print({k:v for k,v in report.items() if k!='checks'},flush=True)
