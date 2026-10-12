"""Literal full original E29 pure cubic, with the true base derivatives."""
from sage.all import *
from pathlib import Path
import time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+750*1024**2,)*2)
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));ed=load(str(p/'independent_ratio_even_derivative_cyclic_direct.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
K=g['coefficient_field'];a,b,h,z=K.gens();P=PolynomialRing(K,names=('l',));l=P.gen();V=PolynomialRing(K,names=tuple('p'+str(i) for i in range(6))+tuple('R'+str(i) for i in range(3)));vv=V.gens()
def cp(f):return V({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
def values(f):
    ans={}
    for ee,c in f.dict().items():
        e=tuple(map(int,ee));j=e[1]+e[4];d=sum(e)
        if e[2] or e[3] or e[7] or j>3 or d>5:continue
        ans[(d,j)]=ans.get((d,j),P.zero())+c*(-1)**e[1]*(-2)**(e[6]+e[8])*l**(e[5]+e[6])
    return ans
def product(A,B):
    return [sum((f*B.get((5-d,j-k),P.zero()) for (d,k),f in A.items()),P.zero()) for j in range(4)]
it=ed['items'][1];i=it['direction'];f=cp(oq['quartic_source_rows'][it['Q_source_index']]);result=vector(P,[0]*4)
for j in range(9):result+=vector(P,product(values(cp(g['first_gauged_fields'][i,j])),values(f.derivative(vv[j]))))
for j,bp in enumerate((h,z)):
    df=V({tuple(map(int,e)):K(c).derivative(bp) for e,c in f.dict().items()})
    result+=vector(P,product(values(cp(g['base_fields'][i,j])),values(df)))
for co,(_,pos,sgn) in zip(it['global_support_coefficient_dicts'],ed['supports'][i]):
    cc=V({tuple(map(int,e)):K(c) for e,c in co.items()})*g['source_unit_coefficients'][i]
    result-=sgn*vector(P,product(values(cc),values(cp(g['six_full_sources'][pos]))))
assert result[0]==0 and result[1]==0 and result[3]==0
save({'ring':P,'full_pure_coefficients':list(result),'record':it,'source_index':29,'scope':__doc__},str(p/'ordinary_pair_plus_zero_E29_cubic.sobj'),compress=False)
print('PASS full E29 pure cubic zero, all material/base; seconds',time.monotonic()-st,flush=True)
