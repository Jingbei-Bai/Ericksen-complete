"""Retain the original opposite offdiagonals on two-pair boundary charts."""
from pathlib import Path
import json,time,sys
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'general_curve_homogeneous_normalized.sobj'));R=d['ring'];K=R.base_ring();H,Z,S,N=K.gens()
lead=load(str(out/'general_curve_dimensionless_leading.sobj'))
T=PolynomialRing(K,names=('x1','y1','x2','y2','x3','y3'));x1,y1,x2,y2,x3,y3=T.gens()
im=[x1,-2*(H*y1+Z*x1)/(1-S),x2,-2*H*(y2+Z*x2)/(1+H*S-Z),x3,2*Z*(y3-H*x3)/((1+H)*S-Z)]
ph=R.hom(im,T);F=[ph(f) for f in d['cubics']];Q=[ph(f) for f in d['quartics']]
def cc(c):return K(c(0,0,H,Z,S,N))
layers=[];records=[]
for missing in range(3):
 active=[j for j in range(3) if j!=missing]
 A=PolynomialRing(K,names=('x','y','z','v'));x,y,z,v=A.gens()
 ims=[A.zero()]*6
 for j,g in enumerate(A.gens()):ims[2*active[j//2]+j%2]=g
 phi=T.hom(ims,A);fs=[phi(f) for f in F if phi(f)];qs=[phi(f) for f in Q]
 B=matrix(K,2,2,[cc(c) for c in lead['blocks'][tuple(active)].list()])
 mix=phi(sum(im[2*active[0]+i]*B[i,j]*im[2*active[1]+j] for i in range(2) for j in range(2)))
 allfs=[mix]+fs+qs
 layers.append({'ring':A,'active':active,'missing':missing,'mixed':mix,'cubics':fs,'quartics':qs,'polynomials':allfs})
 records.append({'missing':missing,'terms':[int(p.number_of_terms()) for p in allfs],'pure_pair_quartic_factors':[[str(A({m:c for m,c in p.dict().items() if (m[2]+m[3]==0 if pair==0 else m[0]+m[1]==0)}).factor()) for pair in range(2)] for p in qs]})
save({'ring':T,'images':im,'cubics':F,'quartics':Q,'layers':layers,'records':records},str(out/'general_curve_two_pair_offdiagonal.sobj'))
(out/'general_curve_two_pair_offdiagonal.json').write_text(json.dumps({'records':records,'elapsed':time.monotonic()-start},indent=2)+'\n')
print(json.dumps({'records':records,'elapsed':time.monotonic()-start}),flush=True)
