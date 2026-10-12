"""Exact two-chart univariate obstructions in the tau normal block.

u=0 is retained separately. On u!=0, the E pure-tau row supplies a
quadratic in u^2,v; its u^2 coefficient is a base-only chart factor.
The resultants below are source-derived, not assertions of global rank.
"""
from sage.all import *
from pathlib import Path
import time,os,resource
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_ratio_p12_axis_tau_elimination.sobj'));U=d['ring'];K=d['coefficient_field'];u,v=U.gens()
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+384*1024**2,)*2)
V=PolynomialRing(K,'tau_v');w=V.gen()
def piece(f,k):return V({int(e[1]):c for e,c in f.dict().items() if int(e[0])==k})
Q=d['E_relations'][0]['after_u'];qu=K(piece(Q,2));qv=piece(Q,0);assert qu and Q==qu*u**2+U(qv(v))
F=d['RA_tau'];G=d['H2_tau'];assert all(e[0] in (0,2,4) for e in F.dict()) and all(e[0] in(0,2) for e in G.dict())
Fp=piece(F,4)*qv**2-piece(F,2)*qv*qu+piece(F,0)*qu**2
Gp=piece(G,0)*qu-piece(G,2)*qv
F0=piece(F,0);G0=piece(G,0)
print('REDUCED degrees',[f.degree() for f in(Fp,Gp,F0,G0)],flush=True)
def sylvester(f,g):
 n,m=f.degree(),g.degree();cols=[]
 for j in range(m):cols.append([f[k-j] if k>=j else K.zero() for k in range(n+m)])
 for j in range(n):cols.append([g[k-j] if k>=j else K.zero() for k in range(n+m)])
 return matrix(K,cols).transpose()
Ks=PolynomialRing(QQ,names=('ta','tr')).fraction_field();a,r=Ks.gens();ph=K.ring().hom([a,a*(1+r)/(1+a*r),1/r,-(1+r)/r],Ks)
def cyc(c):return ph(c.numerator())/ph(c.denominator())
records=[]
for name,ff,gg in [('u_nonzero',Fp,Gp),('u_zero',F0,G0)]:
 if not gg or not ff:
  records.append({'chart':name,'polynomials':(ff,gg),'status':'zero-polynomial; needs another row'});print(name,'zero polynomial',flush=True);continue
 A=sylvester(ff,gg);Ac=A.apply_map(cyc);det=Ac.det();records.append({'chart':name,'polynomials':(ff,gg),'sylvester':A,'cyclic_determinant':det});print(name,'cyc_det_nonzero',bool(det),'shape',A.dimensions(),flush=True)
save({'ring':V,'coefficient_field':K,'u2_coefficient':qu,'u2_free_part':qv,'source_file':'independent_ratio_p12_axis_tau_elimination.sobj','records':records,'scope':__doc__},str(p/'independent_ratio_axis_tau_univariate.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
