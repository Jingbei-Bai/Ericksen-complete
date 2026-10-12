"""Whole-u quadratic remainder after the exact minus-wall range shear.
All original 37 rows and three full pair rows are retained.  This tests
the missing relative u*t^2 / u^2*t^2 generators, rather than freezing u.
"""
from sage.all import *
from pathlib import Path
import json,time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
nr=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));q2=load(str(p/'independent_incoming_simple_wall_quadratics.sobj'));pairs=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+600*1024**2,)*2)
K=pairs['coefficient_field'];a,r,H,Z=K.gens();k=a*r;m=1+k;P=PolynomialRing(K,names=('l','u'));l,u=P.gens();O=nr['ring'];OK=O.base_ring();ph=OK.hom([a,a*(1+r)/m,H/r,(1+r)*Z/r],K)
def mp(ff):return sum((ph(c)*l**int(e[0])*u**int(e[1])*(-2)**int(e[2]) for e,c in O(ff).dict().items()),P.zero())
def rem(ff):return sum((c*(-k/m)**(int(e[0])//2)*l**(int(e[0])%2)*u**int(e[1]) for e,c in P(ff).dict().items()),P.zero())
vals=q2['items'][1]['values'];items=[]
for i in range(37):
    f=mp(vals[i]['second'])-m*mp(nr['matrix'][i,1])-m*l*mp(nr['matrix'][i,2]);rr=rem(f)
    items.append({'label':('old',i),'whole':f,'wall_remainder':rr})
for it in pairs['items']:
    f=it['full_homogeneous_pullback'];out=P.zero()
    for ee,c in f.dict().items():
        e=tuple(map(int,ee));state=e[3:]
        if e[7] or e[8]:continue
        if e[4]==0 and e[5]==0 and e[3]+e[6]==2:weight=(-1)**e[3]
        elif e[4]==1 and sum(state)==1:weight=-m
        elif e[5]==1 and sum(state)==1:weight=-m*l
        else:continue
        out+=weight*K(c)*l**e[0]*u**e[1]*(-2)**e[2]
    items.append({'label':('pair',it['pair']),'whole':out,'wall_remainder':rem(out)})
for it in items:
    if it['wall_remainder']:print(it['label'],'udegrees',sorted({int(e[1]) for e in it['wall_remainder'].dict()}),'value',it['wall_remainder'].factor(),flush=True)
save({'ring':P,'items':items,'scope':__doc__},str(p/'ordinary_pair_minus_whole_u_square.sobj'),compress=False)
print('DONE',sum(bool(it['wall_remainder']) for it in items),'nonzero, seconds',time.monotonic()-st,flush=True)
