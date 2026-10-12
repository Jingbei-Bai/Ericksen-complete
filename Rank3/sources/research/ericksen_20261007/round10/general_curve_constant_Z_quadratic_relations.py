"""Exact cubic-symbol kernel and full quadratic consequences at S=N=0."""
from pathlib import Path
import time,json,sys
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic()
stem='general_curve_proportional_Z' if '--proportional' in sys.argv else 'general_curve_constant_Z'
d=load(str(out/(stem+'_full_sources.sobj')));R=d['ring'];K=R.base_ring();xs=R.gens()
polys=[];names=[]
for pos,pair in enumerate(((0,1),(0,2),(1,2))):
 k=3-sum(pair)
 for j in (2*k,2*k+1):polys.append(d['rows'][pos]*xs[j]);names.append(('mixed',pos,j))
polys.extend(d['rows'][3:7]);names.extend([('cubic',i) for i in range(3,7)])
mons=[xs[i]*xs[2+j]*xs[4+k] for i in range(2) for j in range(2) for k in range(2)]
M=matrix(K,[[p.monomial_coefficient(mon) for p in polys] for mon in mons])
print('matrix built',time.monotonic()-start,flush=True)
ker=M.right_kernel_matrix();print('kernel shape',ker.nrows(),ker.ncols(),time.monotonic()-start,flush=True)
rows=[];denominators=[];checks={}
for i,v in enumerate(ker.rows()):
 den=lcm([c.denominator()/c.denominator().leading_coefficient() for c in v]);den/=den.leading_coefficient();cs=[K(den)*c for c in v]
 f=sum((c*p for c,p in zip(cs,polys)),R.zero());assert f.total_degree()<=2
 rows.append(f);denominators.append(den)
 checks['source_'+str(i)]=not M*vector(cs)
 print('quadratic',i,'terms',f.number_of_terms(),'kernel_denom',str(den.factor()),time.monotonic()-start,flush=True)
save({'ring':R,'tensor_matrix':M,'kernel':ker,'quadratics':rows,'kernel_denominators':denominators,'source_polys':polys,'source_names':names,'checks':checks,'scope':'Exact denominator-cleared quadratic necessary relations on the specified entire straight ratio family. Rank-drop of generic kernel chart must be retained.'},str(out/(stem+'_quadratic_relations.sobj')))
(out/(stem+'_quadratic_relations.json')).write_text(json.dumps({'kernel_dimension':ker.nrows(),'quadratic_terms':[int(f.number_of_terms()) for f in rows],'kernel_denominators':[str(x.factor()) for x in denominators],'checks':checks,'elapsed':time.monotonic()-start},indent=2)+'\n')
