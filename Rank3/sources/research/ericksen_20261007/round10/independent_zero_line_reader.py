"""Original A/Q -> six axis graphs -> contact -> portable source replay."""
from pathlib import Path
from itertools import product
import json,time
from sage.all import load,save,QQ,PolynomialRing,matrix,prod
p=Path(__file__).resolve().parent;start=time.monotonic();checks={}
amp=load(str(p/'general_curve_full_amplitude_forms.sobj'))
raw=load(str(p/'general_curve_two_pair_offdiagonal.sobj'))
axis=load(str(p/'general_curve_rankone_zero_line_coefficient.sobj'))
graph=load(str(p/'general_curve_zero_line_contact_sources.sobj'))
der=load(str(p/'general_curve_zero_line_contact_closure.sobj'))
lift=load(str(p/'general_curve_zero_line_guard_lifts.sobj'))
portable=json.loads((p/'general_curve_zero_line_guard_integer.json').read_text())
R=amp['ring'];K=R.base_ring();C=K.ring();H,Z,S,N=K.gens();xx=R.gens()
old=[C(v) for v in [H,Z,H-Z,Z-1,H+1,H-Z+1,S-1,H*S-Z+1,H*S-Z+S]]
guard=prod(old);A=amp['A'];Q=amp['Q'];odd=amp['odd_pairs']
def unit(v):
    if not v:return False
    return all(C(f).is_constant() or not guard.quo_rem(C(f))[1]
               for a in [K(v).numerator(),K(v).denominator()] for f,e in C(a).factor())
for i in range(6):
    rebuilt=[R.zero()]*3
    for mon,c in raw['cubics'][i].dict().items():
        deg=[sum(mon[2*j:2*j+2]) for j in range(3)];deg[odd[i]]-=1
        assert sorted(deg)==[0,0,2]
        rebuilt[deg.index(2)]+=R({tuple(mon):c})
    checks[f'original_cubic_columns_{i}']=rebuilt==list(A.row(i))
for i in range(3):
    checks[f'original_quartic_reconstruction_{i}']=sum(Q.row(i))==raw['quartics'][i]
    u=xx[2*i]*xx[2*i+1];v=(xx[2*i]+xx[2*i+1])**2
    for j,sgn in enumerate([-QQ(1)/2,-QQ(1)/2,QQ(1)/2]):
        checks[f'pure_quartic_identity_{i}_{j}']=sgn*Q[j,i]==u*v-(2*u**2 if i==j else 0)
checks['six_axis_cases']=[(it['zero_line_index'],it['axis']) for it in axis['items']]==list(product(range(3),range(2)))
for k,it in enumerate(axis['items']):
    i=it['zero_line_index'];side=it['axis'];im=list(xx);im[2*i]=R(side);im[2*i+1]=R(1-side)
    ph=R.hom(im,R);col=[ph(A[j,i]) for j in range(6)]
    checks[f'original_axis_column_{k}']=col==it['column']
    checks[f'axis_odd_rows_zero_{k}']=all(not col[j] for j in range(6) if odd[j]==i)
    dd=[]
    for j in range(3):
        if j==i:continue
        ids=[l for l in range(6) if odd[l]==j]
        B=matrix(K,[[col[l].monomial_coefficient(xx[2*j]),col[l].monomial_coefficient(xx[2*j+1])] for l in ids])
        checks[f'original_axis_block_{k}_{j}']=B==it['blocks'][j] and all(col[l]==B[l0,0]*xx[2*j]+B[l0,1]*xx[2*j+1] for l0,l in enumerate(ids))
        checks[f'original_axis_determinant_{k}_{j}']=B.det()==it['determinants'][j]
        dd.append(B.det())
    g=graph['items'][k];pick=g['picked_det'];nv=K(g['N'])
    alpha=C(dd[pick].numerator()).derivative(C(N))
    checks[f'graph_original_pivot_{k}']=alpha==g['pivot'] and unit(alpha) and unit(dd[pick].denominator())
    checks[f'graph_N_identity_{k}']=not dd[pick].subs({N:nv}) and not nv.derivative(N) and unit(nv.denominator())
    residual=K(dd[1-pick].subs({N:nv}))
    checks[f'other_original_determinant_{k}']=residual==g['other_residual'] and unit(residual/g['primitive_residual'])
    fact=prod((f**e for f,e in g['factors']),C.one())
    checks[f'complete_residual_factorization_{k}']=unit(K(g['primitive_residual'])/fact)
    remaining=[(j,f) for j,(f,e) in enumerate(g['factors']) if not (K(f)==S or K(f)==H*S-Z)]
    checks[f'only_one_nonlinear_factor_{k}']=len(remaining)==1
    jj,f=remaining[0];dl=next(v for v in der['results'] if v['case']==k and v['factor_index']==jj)
    checks[f'contact_initial_factor_{k}']=dl['sources'][0]==f and dl['N']==nv
    for j in range(2):
        v=dl['sources'][j]
        dv=K(v.derivative(C(H)))+S*v.derivative(C(Z))+nv*v.derivative(C(S))
        checks[f'contact_full_chain_rule_{k}_{j}']=unit(dv/dl['sources'][j+1]) and not dl['sources'][j+1].degree(C(N))
    cert=next(v for v in lift['results'] if v['case']==k)
    FR=lift['ring'];hh=FR.base_ring().gen();zz,ss=FR.gens()
    expected=[FR(v(hh,zz,ss,0)) for v in dl['sources'][:3]]
    checks[f'lift_original_contact_sources_{k}']=cert['lift_sources']==expected and len(expected)==3
    checks[f'lift_rational_identity_{k}']=sum((u*v for u,v in zip(cert['cofactors'],expected)),FR.zero())==cert['target']
    native=next(v for v in portable['cases'] if v['case']==k)
    def decode(terms):return C({(int(m[0]),int(m[1]),int(m[2]),0):QQ(c) for m,c in terms})
    checks[f'portable_original_contact_sources_{k}']=[decode(v) for v in native['sources']]==dl['sources'][:3]
    target=decode(native['target_guard']);D=decode(native['denominator_H'])
    checks[f'portable_H_only_denominator_{k}']=bool(D) and all(not D.degree(C(v)) for v in [Z,S,N])
    checks[f'portable_guard_target_{k}']=bool(target) and unit(target) and FR(target(hh,zz,ss,0))==cert['target']
    checks[f'portable_integer_identity_{k}']=sum((decode(u)*decode(v) for u,v in zip(native['cofactors'],native['sources'])),C.zero())==D*target
    print('zero line original case',k,'PASS',round(time.monotonic()-start,2),flush=True)
assert all(checks.values()),[k for k,v in checks.items() if not v]
report={'status':'EXACT_SIX_ZERO_LINE_CASES_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(p/'independent_zero_line_reader.json').write_text(json.dumps(report,indent=2)+'\n')
save({'checks':checks,'scope':'All six indexed-axis cases, original amplitude and conic sources, protected N graph, full actual contact chain rule, and original-source explicit QQ/integer identities. H-only exceptional values are finite and cannot carry an actual interval.'},str(p/'independent_zero_line_reader.sobj'))
print(report['status'],report['check_count'],round(report['elapsed_seconds'],2),flush=True)
