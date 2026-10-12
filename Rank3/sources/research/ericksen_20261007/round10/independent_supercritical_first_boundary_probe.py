"""Exact all-material first supercritical boundary sources, without no-lift claims."""
from sage.all import *
from pathlib import Path
from itertools import combinations
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));K=src['ring'].base_ring();a,r,H,Z=K.gens();q=1+r;alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
R=PolynomialRing(K,names=tuple('sc_p%d'%i for i in range(6))+tuple('sc_v%d'%i for i in range(3)));xs=R.gens()
def cp(f):return R({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
def bi(f,i,j):return R({tuple(map(int,e)):c for e,c in f.dict().items() if sum(e[:6])==i and sum(e[6:])==j})
F12=[bi(cp(f),1,2) for f in src['full_R_sources']]
L=matrix(R,3,3,[cp(f) for f in src['L'].list()])
# First retain every two-RR support with the active-D entries normalized.
TK=PolynomialRing(K,names=('sc_s',)).fraction_field();s=TK.gen();P=PolynomialRing(TK,names=('sc_u','sc_w'));u,w=P.gens()
charts=[([u,w,P(beta/s),0,1,0],[0,P(s),1]),
        ([P(alpha/s),0,u,w,0,1],[P(s),0,1]),
        ([0,P(alpha/(beta*s)),0,1,u,w],[P(s),1,0])]
twos=[]
for j,(pim,vim) in enumerate(charts):
    hom=R.hom(list(map(P,pim+vim)),P);rows=[hom(f) for f in F12]
    ck('two-RR exact L incidence '+str(j),all(hom(f)==0 for f in L*vector(R,xs[6:])))
    constants=[(i,f.constant_coefficient()) for i,f in enumerate(rows) if f and f.total_degree()==0]
    pairs=[]
    for (i,f),(k,g) in combinations(constants,2):
        ff=f.numerator();gg=g.numerator();dd=ff.gcd(gg)
        if dd.degree()==0:
            zz,ca,cb=ff.xgcd(gg);ck('two-RR constant-source Bezout '+str((j,i,k)),ca*ff+cb*gg==zz)
            pairs.append({'rows':[i,k],'sources':[f,g],'polynomial_numerators':[ff,gg],'gcd':zz,'cofactors':[ca,cb]});break
    twos.append({'missing_RR':j,'p_images':pim,'RR_images':vim,'sources':rows,'constant_rows':constants,'constant_coprime_pair':pairs})
    print('TWO SUPPORT',j,'constants',[(i,str(f.factor())) for i,f in constants],'coprime',bool(pairs),flush=True)
# Full RR-torus candidate, retaining every amplitude and all material parameters.
S=PolynomialRing(K,names=('sc_A','sc_B','sc_C','sc_c1','sc_c2'));A,B,C,c1,c2=S.gens()
images=[A,B,C,-c1*c2*B/q,c2*C,c1*A,alpha*c1,beta*c2,S.one()]
ph=R.hom(images,S)
def reduce_roots(f):
    dd={}
    for ex,co in f.dict().items():
        e=list(map(int,ex));n1,n2=e[3]//2,e[4]//2;e[3]%=2;e[4]%=2;co*=(-q/r)**n1*(-q)**n2;ee=tuple(e);dd[ee]=dd.get(ee,K.zero())+co
    return S({e:c for e,c in dd.items() if c})
for j,f in enumerate(L*vector(R,xs[6:])):ck('canonical full torus L '+str(j),reduce_roots(ph(f))==0)
for j,f in enumerate(F12):ck('canonical full torus F12 '+str(j),reduce_roots(ph(f))==0)
h=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));pool=h['polynomial_pool'];old=h['coefficient_field'];km=old.ring().hom([a,a*q/(1+a*r),H/r,q*Z/r],K);cache={}
def cv(co):
    co=old(co)
    if co not in cache:cache[co]=km(co.numerator())/km(co.denominator())
    return cache[co]
H13=[];pulls=[]
for n,rec in enumerate(h['records']):
    dd={}
    for sign,li,ri in rec['full_signed_product_DAG']:
        for ex,c in pool[li].items():
            for ey,b in pool[ri].items():
                ee=tuple(int(x)+int(y) for x,y in zip(ex,ey))
                if sum(ee[:6])!=1 or sum(ee[6:])!=3:continue
                dd[ee]=dd.get(ee,old.zero())+sign*c*b
    f=R({ee:cv(c) for ee,c in dd.items() if c});out=reduce_roots(ph(f));H13.append(f);pulls.append(out)
    print('GENUINE H13',n,'terms',len(f.dict()),'restricted_terms',len(out.dict()),'seconds',round(time.monotonic()-st,2),flush=True)
record={'checks':len(checks),'seconds':time.monotonic()-st,'F12_canonical_all_zero':True,'H13_canonical_nonzero':[j for j,f in enumerate(pulls) if f],'two_RR_constant_coprime_pairs':[bool(it['constant_coprime_pair']) for it in twos],
 'scope':'Exact necessary first supercritical boundary sources. If H13 all zero, records a whole complex all-material candidate family, not an actual or formal lift and not a no-lift result.'}
save({'record':record,'checks':checks,'ring':R,'F12':F12,'H13':H13,'two_RR_supports':twos,'canonical_ring':S,'canonical_images':images,'root_relations':[c1*c1+q/r,c2*c2+q],'H13_canonical_pulls':pulls},str(p/'independent_supercritical_first_boundary_probe.sobj'),compress=False)
(p/'independent_supercritical_first_boundary_probe.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)
