"""Whole incoming source divisibility by l*(v+2) in the axis kernel.
This refines the axis ideal; it is not a frozen v=-2 rank assertion.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent; st=time.monotonic(); checks=[]
ax=load(str(p/'ordinary_single_axis_J9.sobj'))
nr=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
gp=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
U=ax['ring']; K=U.base_ring(); au,av=U.gens(); O=nr['ring']; l,u,v=O.gens()
def ck(n,z):
    assert bool(z),n
    checks.append(n)
def lift(f):return O({(0,int(e[0]),int(e[1])):c for e,c in U(f).dict().items()})
beta=lift(ax['beta']); lam=lift(ax['lambda'])
ck('beta minus two is exactly minus one',beta.subs({v:-2})==-1)
coeff=[]
for i in range(37):
    f=beta*nr['matrix'][i,0]+nr['matrix'][i,3]
    ck('whole l divisibility '+str(i),not f.subs({l:0}))
    ck('whole w divisibility '+str(i),not f.subs({v:-2}))
    z=O(f/(l*(v+2))) if f else O.zero()
    ck('whole exact multiplyback '+str(i),f==l*(v+2)*z)
    coeff.append(z)
vals=[O.one(),O.zero(),O.zero(),O.zero(),O.zero(),l,lam*l,u,v]
def eval_top(ff,degree):
    out=O.zero()
    for ee,c in ff.dict().items():
        e=tuple(map(int,ee))
        if sum(e)!=degree or any(e[j] and not vals[j] for j in range(9)):continue
        out+=K(c)*prod(vals[j]**e[j] for j in range(9) if e[j])
    return out
for j,it in enumerate(gp['items']):
    f=it['full']; xx=f.parent().gens()
    z=beta*eval_top(f.derivative(xx[1]),3)+eval_top(f.derivative(xx[4]),3)
    ck('whole pair l divisibility '+str(j),not z.subs({l:0}))
    ck('whole pair w divisibility '+str(j),not z.subs({v:-2}))
    zz=O(z/(l*(v+2))) if z else O.zero()
    ck('whole pair exact multiplyback '+str(j),z==l*(v+2)*zz)
    coeff.append(zz)
def ataxis(f):return U({(int(e[1]),int(e[2])):c for e,c in O(f).dict().items() if not e[0]})
B=matrix(U,ax['matrix'])
for i,z in enumerate(coeff):
    bz=ataxis(z)
    ck('refined l w k axis coefficient '+str(i),B[i,5]==(av+2)*bz)
    B[i,5]=bz
C=B.apply_map(lambda f:K(f.subs({au:0,av:-2})))
KB=FractionField(PolynomialRing(QQ,names=('a','r')));a,r=KB.gens();q=1+r;m=1+a*r
ph=K.hom([a,a*q/m,1/r,-q/r],KB)
D=C.apply_map(ph)
rank=D.rank(); piv=D.transpose().pivots()
print('CYCLIC REFINED RANK',rank,'PIVOTS',piv,flush=True)
for i,row in enumerate(D.rows()):
    if any(row):print('SUP',i,[j for j,z in enumerate(row) if z],flush=True)
dd=D.matrix_from_rows(piv).det() if rank==9 else KB.zero()
print('TRUE NINE MINOR',dd.factor() if dd else 0,flush=True)
save({'ring':U,'whole_ring':O,'whole_divided_coefficients':coeff,'matrix':B,'minus_u_zero_matrix':C,'cyclic_matrix':D,'cyclic_rank':rank,'pivot_rows':piv,'cyclic_minor':dd,'checks':checks,'columns':('range0','n21','n23','tau','u_eta','l_w_k','eta2','eta_k','k2')},str(p/'ordinary_axis_v_minus_relative_factor.sobj'),compress=False)
rep={'status':'PASS WHOLE L*(V+2) DIVISIBILITY AND COMPLETE MATRIX','checks':len(checks),'cyclic_rank':rank,'minor':str(dd.factor()) if dd else '0','seconds':time.monotonic()-st}
(p/'ordinary_axis_v_minus_relative_factor.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep,flush=True)
