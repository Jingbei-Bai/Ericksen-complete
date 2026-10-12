"""Whole 28+9 source second obstruction in all conic endpoint tangents."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_full_regular_chart_H2.sobj'));g=load(str(p/'independent_subcritical_single_RR_full_complete_E_chart_H2.sobj'));P=d['ring'];OK=P.base_ring();kk,rr=OK.gens();nn=load(str(p/'independent_subcritical_single_RR_conic_endpoint_module.sobj'));E=nn['ring'];K=E.base_ring();k,r=K.gens();xb=E.gen();ph=OK.ring().hom([k,r],K);T=PolynomialRing(E,names=('qa_v','qa_s','qa_w'));v,s,w=T.gens();x=T(xb)-v/24-s/2;z=-x+s;zero=T.zero()
def mp(f,imgs):return sum((E(ph(c.numerator())/ph(c.denominator()))*prod(iv**int(n) for iv,n in zip(imgs,ee)) for ee,c in f.dict().items()),T.zero())
bd=load(str(p/'independent_critical_single_RR3_normal.sobj'));BP=bd['ring'];BK=BP.base_ring();bh=BK.ring().hom([k/r,r,2,-1],K)
def bm(f):return sum((E(bh(c.numerator())/bh(c.denominator()))*prod(iv**int(n) for iv,n in zip((x,z,w,zero,zero,zero,zero),ee)) for ee,c in f.dict().items()),T.zero())
graphs=list(map(bm,bd['graph_images']));images=(x,z,w,graphs[3],graphs[4],graphs[5]+v,graphs[6],zero,zero);sources=d['sources']+[it['source'] for it in g['items']];polys=[mp(f,images) for f in sources];base=(E(xb),E(-xb),E.zero(),E.zero(),E.zero(),E(2*xb),E.zero(),E.zero(),E.zero())
def em(f):return sum((E(ph(c.numerator())/ph(c.denominator()))*prod(iv**int(n) for iv,n in zip(base,ee)) for ee,c in f.dict().items() if not any(n and not iv for n,iv in zip(ee,base))),E.zero())
M=matrix(E,len(sources),9,lambda i,j:em(sources[i].derivative(P.gen(j))));print('M rank',M.rank(),'sec',time.monotonic()-st,flush=True)
for i,f in enumerate(polys):
 for mv in (T.one(),v,s,w):assert f.monomial_coefficient(mv)==0,(i,mv)
mons=(v*v,v*s,s*s,v*w,s*w,w*w);Q=matrix(E,len(sources),6,lambda i,j:polys[i].monomial_coefficient(mons[j]));left=M.left_kernel().basis_matrix();obs=left*Q;nonzero=[row for row in obs if any(row)];print('obs rank',obs.rank(),'nonzero',len(nonzero),'sec',time.monotonic()-st,flush=True)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 n,d=c.numerator(),c.denominator();gg=n.gcd(d);n,d=n//gg,d//gg;lc=d.leading_coefficient();return K(n/lc)/K(d/lc)
def ne(c):return E([nk(u) for u in E(c).list()])
clean=[vector(E,[ne(c) for c in row]) for row in nonzero]
for row in clean:print('obstruction',sum((c*mv for c,mv in zip(row,mons)),T.zero()),flush=True)
save({'ring':T,'sources':sources,'normal_matrix':M,'images':images,'boundary_pullbacks':polys,'quadratic_matrix':Q,'left_kernel':left,'quadratic_obstructions':clean,'scope':__doc__},str(p/'independent_subcritical_single_RR_conic_endpoint_complete_quadratic.sobj'),compress=False)
