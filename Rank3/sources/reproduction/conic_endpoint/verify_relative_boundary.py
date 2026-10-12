"""Reconstruct complete kappa=t=0 boundary directly from original sources.

No saved boundary/graph matrix is read. Parameters k,r,H,Z remain symbolic.
All identities are polynomial identities before local completion.
"""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent
p=out/'source_inputs'
start=time.monotonic()
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
old=full['coefficient_field']
OK=qf['ring'].base_ring()
K=PolynomialRing(QQ,names=('k','r','H','Z')).fraction_field()
k,r,H,Z=K.gens(); m=1+k;q=1+r
P=PolynomialRing(K,names=('x','z','w'));x,z,w=P.gens();s=x+z
images=(x,P.one(),z,w,P.zero(),P.zero(),P.zero(),P.zero(),P.one())
ph=old.ring().hom([k/r,k*q/(r*m),H/r,q*Z/r],K)
po=OK.ring().hom([k/r,r,H,Z],K)
cache={}; denoms=set(); checks=[]
def cv(c,F):
    key=(str(F),c)
    if key not in cache:
        hh=ph if F==old else po
        cc=hh(c.numerator())/hh(c.denominator())
        cache[key]=cc;denoms.add(cc.denominator())
    return cache[key]
def mpdict(dd,F,layer=None):
    ans=P.zero()
    for ee,c in dd.items():
        ee=tuple(map(int,ee))+(0,)*(9-len(ee))
        if layer is not None and (sum(ee[:6]),sum(ee[6:]))!=layer: continue
        if any(n and not image for n,image in zip(ee,images)):continue
        ans+=cv(c,F)*prod(image**n for n,image in zip(ee,images) if n)
    return ans
def mp(f,layer=None):return mpdict(f.dict(),f.parent().base_ring(),layer)
def ck(label,test):
    assert bool(test),label
    checks.append(label)
f0=((2*H*Z-2*H-4*Z**2+2*Z)*x+(4*H*Z-2*H-4*Z**2+2*Z)*z)*w/(H*Z)
f1=-2*Z*(H-1)*w/(H*(Z-1))
f2=((2*H*Z-2*H-4*Z**2+4*Z)*x+(2*H*Z-2*H-4*Z**2+2*Z)*z)/Z
d=4*r*(H+1)**3*s*w/(q*H**2)
es=[f0,f1,f2]
DC=matrix(qf['D'])*matrix(qf['C'])
graph=[]
for i,f in enumerate(qf['full_R_sources']):
    g=mp(f,(4,0))+d*mp(f,(1,2))+sum(mp(DC[i,j])*es[j] for j in range(3))
    ck('R'+str(i)+' full relative graph identity',g==0);graph.append(g)
g=mp(qf['full_Q0'],(2,1))+sum(mp(qf['Q0_L_cofactors'][j])*es[j] for j in range(3))
ck('Q0 full relative graph identity',g==0);graph.append(g)
pool=full['polynomial_pool']
def product_layer(record,layer):
    ans=P.zero()
    for sign,li,ri in record['full_signed_product_DAG']:
        for le,lc in pool[li].items():
            for re,rc in pool[ri].items():
                ee=tuple(int(a)+int(b) for a,b in zip(le,re))
                if (sum(ee[:6]),sum(ee[6:]))!=layer:continue
                if any(n and not image for n,image in zip(ee,images)):continue
                ans+=sign*cv(lc,old)*cv(rc,old)*prod(image**n for n,image in zip(ee,images) if n)
    return ans
for i,item in enumerate(hs['items']):
    g=product_layer(full['records'][i],(4,1))+d*product_layer(full['records'][i],(1,3))
    g+=sum(mpdict(item['L_cofactors_original'][j],old)*es[j] for j in range(3))
    ck('H'+str(i)+' full relative graph identity',g==0);graph.append(g)
    print('H',i,'PASS',round(time.monotonic()-start,2),flush=True)
for item in ev['items']:
    g=mp(item['pure_p4'])+sum(mp(item['A_p_quadratic'][j])*es[j] for j in range(3))
    g+=cv(item['square_unit'],OK)*es[int(item['square_component'])]**2
    graph.append(g)
C=2*k/m**2*(x**3*z+x*z**3)+(10*x*x-4*x*z-6*z*z)*w*w/m+k*(-6*x*x-4*x*z+10*z*z)/m
ck('Q1 exact s-square division',graph[25]==s*s*(-2*x*z/m-6*w*w/k+26))
ck('Q5 exact s-square division',graph[26]==s*s*(-2*k*x*z/m+26*w*w-6*k))
ck('Q6 exact C graph relation',graph[27]==C)
sub=P.hom([x,-x,w],P)
ck('C conic exact factorization',sub(C)==-4*k*x*x/m**2*(x*x-2*m*(k+w*w)/k))
cx=C.derivative(x)-C.derivative(z)
ratio=(3*w*w-k)/(2*(w*w+k))
test=sub(C.derivative(z))+sub(cx)*ratio
conic=x*x-2*m*(k+w*w)/k
ck('nonreduced conic lift',test.numerator().reduce(P.ideal(conic).groebner_basis())==0)
center_den_factors=set()
for den in denoms:
    dd=den(H=1,Z=-1)
    ck('source coefficient regular at cyclic base '+str(den),bool(dd))
    for factor,power in dd.factor():center_den_factors.add(str(factor))
payload={'status':'PASS','checks':len(checks),'identities':checks,'center_denominator_factors':sorted(center_den_factors),
    'graph':list(map(str,[f0,f1,f2,d])),'C':str(C),'seconds':time.monotonic()-start,
    'scope':__doc__}
(out/'relative_boundary_reconstruction.json').write_text(json.dumps(payload,indent=2)+'\n')
save({'ring':P,'graph':[f0,f1,f2,d],'C':C,'graph_residuals':graph,'checks':checks},str(out/'relative_boundary_reconstruction.sobj'),compress=False)
print('PASS',len(checks),'DENOMINATOR FACTORS',sorted(center_den_factors),flush=True)
