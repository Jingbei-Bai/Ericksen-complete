"""Uniform rank-two certificate for the critical line amplitude equations."""
from sage.all import *
from pathlib import Path
import time,json,itertools
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_line_delta_coefficients.sobj'));K=d['field'];a,r,t=K.gens();q=1+r;M=d['matrix'];w=vector(K,[-r*(t+1)**2,q*(t+1)**2,r*q*t*t]);assert M*w==0;N=M*diagonal_matrix(K,w);T=PolynomialRing(QQ,'tt');tt=T.gen();ph=K.ring().hom([T.zero(),T.zero(),tt],T);rows=[];items=[]
for i in range(N.nrows()):
 v=list(N.row(i));den=lcm([f.denominator() for f in v]);pp=[K.ring()(f*den) for f in v];g=gcd(pp);u=K(g.lc())
 for f,e in g.factor():
  if f.degree(K.ring().gens()[2])==0:
   assert any(K(f)/v0 in QQ and K(f)/v0!=0 for v0 in (a,r,q,1+a*r)),(i,f)
   u*=K(f/f.lc())**e
 out=[K(f)/u for f in pp];assert all(v0.derivative(a)==0 and v0.derivative(r)==0 for v0 in out)
 rr=tuple(ph(v0.numerator())/QQ(v0.denominator()) for v0 in out);items.append({'row':i,'denominator':den,'material_unit':u,'normalized_row':rr})
 if rr not in rows:rows.append(rr)
g=T.zero();co=[];sources=[];links=[]
for i,j in itertools.combinations(range(len(rows)),2):
 f=rows[i][0]*rows[j][1]-rows[i][1]*rows[j][0]
 if not f:continue
 ng,aa,bb=g.xgcd(f);co=[aa*v for v in co]+[bb];sources.append(f);links.append((i,j));g=ng
 # Protected t(t+1) support: once no other factors remain this is complete.
 if all(f==tt or f==tt+1 for f,e in g.factor()):break
assert sum((v*f for v,f in zip(co,sources)),T.zero())==g
print('GCD',g.factor(),'minors',len(sources),'rows',len(rows),'seconds',time.monotonic()-st,flush=True)
save({'field':K,'short_kernel':w,'all_rows_kernel_check':True,'transform_items':items,'ring':T,'distinct_rows':rows,'minor_sources':sources,'minor_links':links,'cofactors':co,'target':g,'scope':__doc__},str(p/'independent_delta_critical_line_rank_certificate.sobj'),compress=False)
(p/'independent_delta_critical_line_rank_certificate.json').write_text(json.dumps({'rows':len(items),'minor_sources':len(sources),'gcd':str(g.factor()),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
