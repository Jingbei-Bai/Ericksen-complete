from sage.all import *
from pathlib import Path
OUT=Path(__file__).absolute().parent;d=load(str(OUT/'homothety_invariant_linear_chart.sobj'));K=d['source']['ring'].base_ring();a,b,s,t=K.gens();BR=K.ring();alpha=b*(a-1);beta=a*(b-1);gamma=a-b
ABC=PolynomialRing(QQ,names=('X','Y','Z'));X,Y,Z=ABC.gens()
def compress(poly,degree):
 out=ABC(0)
 for (ea,eb,es,et),cf in poly.dict().items():
  if es%2 or et%2:raise ValueError('odd')
  # Group a,b coefficients first.
 groups={}
 for mon,cf in poly.dict().items():
  ea,eb,es,et=mon;key=(es//2,et//2);groups[key]=groups.get(key,K(0))+K(cf)*a**ea*b**eb
 for (m,n),co in groups.items():
  i=m+n-degree;j=degree-n;k=degree-m
  if min(i,j,k)<0:raise ValueError((m,n,degree))
  val=co/(alpha**i*beta**j*gamma**k)
  out+=QQ(val)*X**i*Y**j*Z**k
 return out
facs=list(d['quadratic_fiber_leaders'][0].numerator().factor());G=max((BR(ff) for ff,n in facs),key=lambda z:len(z.dict()))
facs=list(d['quadratic_fiber_leaders'][1].numerator().factor());F=max((BR(ff) for ff,n in facs),key=lambda z:len(z.dict()))
GG=compress(G,4);FF=compress(F,2)
print('F=',FF,'G=',GG,'G-F²=',(GG-FF**2).factor(),flush=True)
save({'F':FF,'G':GG,'difference':GG-FF**2},str(OUT/'homothety_fiber_leading_compression.sobj'))
