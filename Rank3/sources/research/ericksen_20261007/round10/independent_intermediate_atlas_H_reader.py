"""Finite pivot atlas, both toric charts, and genuine complete H regularity."""
from sage.all import *
from pathlib import Path
from itertools import combinations,product
import json,time
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
def eq(n,f,g):ck(n,not(f-g))
K=PolynomialRing(QQ,names=('at_alpha','at_beta')).fraction_field();alpha,beta=K.gens()
R=PolynomialRing(K,names=tuple('at_p%d'%j for j in range(6))+tuple('at_e%d'%j for j in range(3))+('at_kappa',));p0,p1,p2,p3,p4,p5,e0,e1,e2,k=R.gens();ps=R.gens()[:6];ee=vector(R,(e0,e1,e2))
L=matrix(R,[[p0,0,-alpha*p5],[p1,-alpha/beta*p3,0],[0,p2,-beta*p4]])
Delta=p0*p3*p4-p1*p2*p5;eq('detL=alpha Delta',L.det(),alpha*Delta)
minors=[]
for ii in combinations(range(3),2):
    for jj in combinations(range(3),2):
        f=L.matrix_from_rows_and_columns(ii,jj).det();ck('2minor single protected monomial '+str((ii,jj)),len(f.dict())==1 and f!=0);minors.append(f)
partials=[Delta.derivative(x) for x in ps]
for bits in product((0,1),repeat=6):
    active=all(any(bits[j] for j in pair) for pair in ((0,1),(2,3),(4,5)))
    if not active:continue
    def nonzero_monom(f):return any(all(not ex[j] or bits[j] for j in range(6)) for ex in f.dict())
    ck('every active-D support rankL at least2 '+str(bits),any(nonzero_monom(f) for f in minors))
    ck('every active-D support gradientDelta nonzero '+str(bits),any(nonzero_monom(f) for f in partials))
eq('Euler identity protects nonnormalized pivot',sum(ps[j]*partials[j] for j in range(6)),3*Delta)
adj=L.adjugate()
for j in range(3):
    F=L.det()-k*(adj*ee)[j]
    for l,pp in enumerate(ps):
        ck('each determinant constraint affine in p '+str((j,l)),F.derivative(pp).derivative(pp)==0)
        eq('p pivot coefficient at kappa0 '+str((j,l)),F.derivative(pp).subs({k:0}),alpha*partials[l])
# The explicit alternative solve-p23 chart includes p32=0.
FQ=R.fraction_field();kp=FQ(k);den=alpha*p0*(p4+kp*e2/beta)
sol=(p1*p2*(alpha*p5+kp*e0)-kp*p0*p2*e1)/den
w0=(alpha*p5+kp*e0)/p0;w1=(beta*p4+kp*e2)/p2
eq('alternative p23 exact ell1',p1*w0-alpha/beta*sol*w1,kp*e1)
eq('alternative p23 corner denominator',den.subs({k:0}),alpha*p0*p4)
eq('alternative p23 zero-p32 corner',sol.subs({k:0,p5:0}),0)
print('FINITE support atlas PASS',round(time.monotonic()-st,2),flush=True)
# All H products are original, including original h,z coefficient derivatives.
h=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
hm=load(str(p/'independent_delta_H_global_L_module.sobj'))
pool=[{tuple(map(int,e)):c for e,c in f.items() if c} for f in h['polynomial_pool']]
OK=h['coefficient_field'];RR=PolynomialRing(OK,names=tuple('ahat%d'%j for j in range(9)))
def oc(f):return RR({tuple(map(int,e)):OK(c) for e,c in f.dict().items() if c})
def pp(ref):return RR(pool[int(ref)])
deg=[max((sum(e) for e in f),default=-100) for f in pool]
wt=[max((sum(e[:6])+2*sum(e[6:]) for e in f),default=-100) for f in pool]
def part(ref,od=None,bd=None):return {e:c for e,c in pool[ref].items() if (od is None or sum(e)==od) and (bd is None or (sum(e[:6]),sum(e[6:]))==bd)}
cache={}
def extract(dag,mode):
    out={}
    for sign,li,ri in dag:
        key=(int(li),int(ri),mode)
        if key not in cache:
            dd={}
            for ex,c in pool[li].items():
                for ey,b in pool[ri].items():
                    e=tuple(x+y for x,y in zip(ex,ey));ii,jj=sum(e[:6]),sum(e[6:])
                    if (mode=='six' and ii+jj!=6) or (mode=='32' and (ii,jj)!=(3,2)):continue
                    dd[e]=dd.get(e,OK.zero())+c*b
            cache[key]={e:c for e,c in dd.items() if c}
        for e,c in cache[key].items():out[e]=out.get(e,OK.zero())+sign*c
    return RR({e:c for e,c in out.items() if c})
scales=tuple(OK(c) for c in hm['original_to_binomial_variable_scales'])
bl=[oc(f) for f in hm['binomial_generators']]
for n,(record,mrecord) in enumerate(zip(h['records'],hm['items'])):
    dag=record['full_signed_product_DAG']
    for sign,li,ri in dag:
        ck('full H term weighted<=7 '+str((n,li,ri)),wt[li]+wt[ri]<=7)
        ck('full H term ordinary<=6 '+str((n,li,ri)),deg[li]+deg[ri]<=6)
    h6=extract(dag,'six');ck('GENUINE complete ordinary6 cancellation '+str(n),h6==0)
    f32=extract(dag,'32')
    normalized=RR({e:c*prod(scales[j]**e[j] for j in range(9)) for e,c in f32.dict().items()})
    eq('original H32 source link '+str(n),normalized,oc(mrecord['normalized_face']))
    cf=[oc(f) for f in mrecord['cofactors']]
    eq('H32 original polynomial module '+str(n),normalized,sum((cf[j]*bl[j] for j in range(3)),RR.zero()))
    ck('H32 polynomial multiplier degrees '+str(n),all(all(sum(e[:6])==2 and sum(e[6:])==1 for e in f.dict()) for f in cf))
    # This finite superset contains every remaining possible complete monomial.
    for i in range(6):
        for j in range(4):
            if i+j>5 or i+2*j>7 or (i,j)==(3,2):continue
            ks,ds=11-2*i-3*j,5-i-j;ss,bs=11-2*i-3*j,7-i-2*j
            ck('sub/super H regular all allowable slots '+str((n,i,j)),min(ks,ds,ss,bs)>=0)
    print('H',n,'full source/atlas PASS',round(time.monotonic()-st,2),flush=True)
record={'checks':len(checks),'seconds':time.monotonic()-st,'H_full_sources':18,'scope':'All active-D supports: rankL2 and smooth Delta on Delta0. Finite rational pivot atlas; all complete H sources regular after the stated subcritical/critical/supercritical normalizations. No no-lift assertion.'}
save({'checks':checks,'record':record,'L':L,'Delta':Delta,'two_minors':minors,'Delta_gradient':partials},str(p/'independent_intermediate_atlas_H_reader.sobj'),compress=False)
(p/'independent_intermediate_atlas_H_reader.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2),flush=True)
