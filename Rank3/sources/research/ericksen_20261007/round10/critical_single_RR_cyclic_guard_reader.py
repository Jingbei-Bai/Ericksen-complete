"""All original critical single-RR support cofactors are regular at cyclic base."""
from sage.all import *
from pathlib import Path
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_single_RR_candidate_certificate.sobj'));OP=d['ring'];OK=OP.base_ring();K=PolynomialRing(QQ,names=('a','r')).fraction_field();a,r=K.gens();ph=OK.ring().hom([a,r,K(1),K(-1)],K);P=PolynomialRing(K,names=OP.variable_names());checks=[]
def ck(n,v):assert bool(v),n;checks.append(n)
def cf(c):return ph(c.numerator())/ph(c.denominator())
def ff(f):return P({tuple(map(int,e)):cf(c) for e,c in f.dict().items()})
for i,record in enumerate(d['denominator_audit']):
 val=ph(record['factor']);ck('cyclic denominator nonzero '+str(i),val!=0)
 for f,e in val.numerator().factor():
  f=K.ring()(f)
  ck('material cyclic denominator protected '+str((i,f)),f.total_degree()==0 or f/f.leading_coefficient() in (a,r,a-1,r+1,1+a*r))
for it in d['linear_delta_lifts']:
 ck('original cyclic support delta identity '+str((it['axis'],it['mask'])),sum((cf(c)*ff(f) for c,f in zip(it['cofactors'],it['rows'])),P.zero())==P.gen(7))
for it in d['source_graph_lifts']:
 for j,(co,target) in enumerate(zip(it['echelon_cofactors'],it['echelon'])):
  ck('original cyclic complete graph identity '+str((it['axis'],it['mask'],j)),sum((cf(c)*ff(f) for c,f in zip(co,it['rows'])),P.zero())==ff(target))
ck('all 45 supports remain covered',{(it['axis'],it['mask']) for it in d['linear_delta_lifts']+d['source_graph_lifts']}=={(i,j) for i in range(3) for j in range(1,16)})
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__};save({'checks':checks,'scope':__doc__},str(p/'critical_single_RR_cyclic_guard_reader.sobj'),compress=False);(p/'critical_single_RR_cyclic_guard_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
