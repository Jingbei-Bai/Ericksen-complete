"""Exact two-column rank for the nilpotent Schur no-lift argument.
The ninth generator eta*k need not be eliminated linearly.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'ordinary_axis_v_minus_relative_factor.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
U=d['ring'];K=U.base_ring();R=PolynomialRing(K,'v');v=R.gen()
def cv(f):return R({int(e[1]):c for e,c in U(f).dict().items() if not e[0]})
M=d['matrix'].apply_map(cv)
def ck(n,z):assert bool(z),n;checks.append(n)
arows=(7,11,25,37);brows=(1,2,3,34,35,36)
for i in arows:ck('A support '+str(i),all(not M[i,j] for j in (0,1,2,4,6,8)))
for i in brows:ck('B support '+str(i),all(not M[i,j] for j in (0,2,3,4,5,7)))
items=[]
for i,j in combinations(arows,2):
    f=M[i,3]*M[j,5]-M[i,5]*M[j,3]
    if f:items.append({'rows':(i,j),'polynomial':f});print('PAIR',i,j,'DEGREE',f.degree(),flush=True)
KK=FractionField(PolynomialRing(QQ,names=('a','r')));a,r=KK.gens();q=1+r;m=1+a*r
ph=K.hom([a,a*q/m,1/r,-q/r],KK);S=PolynomialRing(KK,'v');sv=S.gen()
def cyc(f):return S({int(e):ph(c) for e,c in R(f).dict().items()})
for it in items:
    it['cyclic']=cyc(it['polynomial'])
    print('CYCLIC',it['rows'],it['cyclic'].factor(),flush=True)
gc=gcd([it['cyclic'] for it in items]);print('CYCLIC GCD',gc,flush=True)
save({'ring':R,'matrix':M,'items':items,'cyclic_ring':S,'cyclic_gcd':gc,'checks':checks},str(p/'ordinary_axis_nilpotent_schur_two_minors.sobj'),compress=False)
for i,j in combinations(range(len(items)),2):
    ff,gg=items[i]['cyclic'],items[j]['cyclic']
    if gcd(ff,gg).degree()<=0:
        res=ff.resultant(gg)
        print('COPRIME PAIR',items[i]['rows'],items[j]['rows'],'RESULTANT',res.factor(),flush=True)
        save({'ring':R,'matrix':M,'items':items,'cyclic_ring':S,'cyclic_gcd':gc,'coprime_pair':(i,j),'cyclic_resultant':res,'checks':checks},str(p/'ordinary_axis_nilpotent_schur_two_minors.sobj'),compress=False)
        break
rep={'status':'EXACT NILPOTENT SCHUR TWO-MINORS','checks':len(checks),'cyclic_gcd':str(gc),'seconds':time.monotonic()-st}
(p/'ordinary_axis_nilpotent_schur_two_minors.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep,flush=True)
