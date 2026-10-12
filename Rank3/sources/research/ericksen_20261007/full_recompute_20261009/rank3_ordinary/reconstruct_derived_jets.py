"""Derive axis and extra-tau-wall jets from the freshly rebuilt full jets.

These are exact polynomial changes of the normal variables, with all base
and tangent parameters retained. No historical jet coefficients are read.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations_with_replacement
import json,time
OUT=Path(__file__).absolute().parent;DATA=OUT/'rebuild/round10'
first=load(str(DATA/'independent_ratio_incoming12_generic_normal.sobj'))
second=load(str(DATA/'independent_ratio_incoming12_generic_normal_two.sobj'))
F=first['ring'];K=F.base_ring();a,b,h,z=K.gens();lf,uf,vf=F.gens()
M=first['matrix'];C=second['coefficient_matrix'];nx=second['normal_exponents']
assert M.nrows()==C.nrows()==37
assert first['lambda'].derivative(lf)==0
checks=[];start=time.monotonic()
def ck(name,ok):
    assert ok,name
    checks.append(name)

ck('incoming lambda on null wall',first['lambda'].subs({uf:0,vf:-2})==-2)
def exponents(n):
    result=[]
    for i in range(n):
        e=[0]*n;e[i]=1;result.append(tuple(e))
    for i,j in combinations_with_replacement(range(n),2):
        e=[0]*n;e[i]+=1;e[j]+=1;result.append(tuple(e))
    return result

# Axis: l becomes a seventh normal; RR1 = eta + lambda(v)*l.
U=PolynomialRing(K,names=('relax_u','relax_v'));u,v=U.gens()
def at_axis(f):
    return U({(int(e[1]),int(e[2])):c for e,c in F(f).dict().items() if not e[0]})
V=PolynomialRing(U,names=tuple('axis_n'+str(i) for i in range(7)))
x=V.gens();lam=at_axis(first['lambda'])
images=(x[0],x[1],x[2],x[3],x[5]-lam*x[4],x[6])
ax=exponents(7);rows=[]
for i in range(37):
    value=sum((V(at_axis(M[i,j]))+V(at_axis(M[i,j].derivative(lf)))*x[4])*images[j] for j in range(6))
    for j,e in enumerate(nx):
        if sum(e)==2:value+=V(at_axis(C[i,j]))*prod(images[k]**int(e[k]) for k in range(6))
    ck('axis exact degree at most two '+str(i),value.total_degree()<=2)
    rows.append([value.monomial_coefficient(prod(x[j]**e[j] for j in range(7))) for e in ax])
matrix_axis=matrix(U,rows)
save(dict(coefficient_field=K,ring=U,normal_positions=(1,2,3,4,5,6,'tau'),
          source_records=first['source_records'],normal_exponents=ax,
          coefficient_matrix=matrix_axis,linear_matrix=matrix_axis.matrix_from_columns(range(7)),
          scope='Exact polynomial transport of complete fresh incoming jets, including l derivatives of every linear coefficient.'),
     str(DATA/'independent_ratio_p12_axis_generic_normal_two.sobj'),compress=False)
print('Fresh axis 37 by 35 coefficients reconstructed.',flush=True)

# Extra-tau wall: retain its quadratic amplitude algebra and all tangents.
U=PolynomialRing(K,'gnq_l');l=U.gen()
dd=-(a-b)*z*z/(b*(a-1));mod=l*l-dd
r=(a-b)/(a*(b-1));q=1+r;Z=r*z/q;kap=(Z-1)/Z
def ev(f):
    f=F(F(f).subs({uf:0,vf:-2}))
    return U({int(e[0]):c for e,c in f.dict().items()})%mod
V=PolynomialRing(U,names=('kernel','eta','tau'));t,e,tau=V.gens()
images=(-t+kap*l*tau,V.zero(),V.zero(),t,e,tau)
qx=[ee for ee in exponents(3) if sum(ee)==2]
items=[]
for i in range(37):
    linear=sum(V(ev(M[i,j]))*images[j] for j in range(6))
    ck('extra-tau all three true kernel directions '+str(i),all(c%mod==0 for c in linear.coefficients()))
    quadratic=V.zero()
    for j,ee in enumerate(nx):
        if sum(ee)==2:quadratic+=V(ev(C[i,j]))*prod(images[k]**int(ee[k]) for k in range(6))
    rhs=[quadratic.monomial_coefficient(t**ee[0]*e**ee[1]*tau**ee[2])%mod for ee in qx]
    mr=M.row(i)
    for tangent in (lf,uf,vf):
        rhs.extend([ev((mr[3]-mr[0]).derivative(tangent)),ev(mr[4].derivative(tangent)),
                    (kap*l*ev(mr[0].derivative(tangent))+ev(mr[5].derivative(tangent)))%mod])
    items.append(dict(source_row=i,linear_row=[ev(f) for f in mr],rhs=rhs))
payload=dict(ring=U,relation=mod,kappa=kap,quadratic_exponents=qx,
             tangent_order=['dl','du','dv'],items=items,complete=True,
             scope='All 37 true second jets reconstructed by exact kernel substitution in freshly generated normal jets, retaining dl,du,dv.')
save(payload,str(DATA/'independent_ratio_incoming12_generic_nullwall_all.sobj'),compress=False)
selected={0,1,2,5,17,29,35};payload=dict(payload,items=[it for it in items if it['source_row'] in selected])
save(payload,str(DATA/'independent_ratio_incoming12_generic_nullwall_even.sobj'),compress=False)
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-start,
            fresh_sources_only=True,historical_coefficient_inputs=0)
(OUT/'derived_jet_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS',len(checks),'derived-jet identities',flush=True)
