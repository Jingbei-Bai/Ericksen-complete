"""Verify that the graph used by the cyclic quartic reader solves actual R tops."""
from sage.all import *
from pathlib import Path
import json,hashlib
out=Path(__file__).absolute().parent;research=out.parent.parent;rank=out.parent/'rank3';checks=[];imports=[]
def rd(p):imports.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()));return load(str(p))
def ck(n,v):assert bool(v),n;checks.append(n)
g=rd(research/'round10/independent_ratio_p_graph_original_sources_cyclic.sobj');f=rd(rank/'fresh_global_polynomials.sobj');T=g['p_ring'];K=T.base_ring();a,r=K.gens();q=1+r;m=1+a*r;OK=f['ring'].base_ring();ph=OK.ring().hom([a,a*q/m,1/r,-q/r],K);cache={}
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
den=g['graph_denominator'];ns=g['graph_numerators'];ck('graph denominator nonzero',den!=0)
def pull(f):
 out=T.zero()
 for e,c in f.dict().items():
  if sum(e)!=4:continue
  re=sum(e[6:]);ck('actual ordinary R top at most linear RR',re<=1)
  out+=T({tuple(map(int,e[:6])):cv(c)})*den**(1-re)*prod(ns[j]**int(e[6+j]) for j in range(3))
 return out
for j,source in enumerate(f['sources']):ck('graph solves fresh complete ordinary R top '+str(j),pull(source)==0)
def pure(p):return T({tuple(map(int,e[:6])):cv(c) for e,c in p.dict().items()})
C=matrix(T,3,3,[pure(p) for p in f['C'].list()]);L=matrix(T,3,3,[pure(p) for p in f['L'].list()]);ck('generic graph RR coefficient invertible',C.det()*L.det()!=0)
report=dict(status='PASS',count=len(checks),checks=checks,imports=imports,scope='The candidate rational graph used in cyclic quartic pullbacks satisfies all six freshly derived ordinary-leading R equations; its RR coefficient is generically invertible, identifying the unique source graph.')
(out/'cyclic_graph_source.json').write_text(json.dumps(report,indent=2));print('PASS',len(checks))
