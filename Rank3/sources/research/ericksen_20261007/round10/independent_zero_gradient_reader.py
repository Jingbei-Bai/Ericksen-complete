"""All seven zero-gradient direction patterns: original-source audit."""
from pathlib import Path
from itertools import product,combinations
import time,json
from sage.all import load,save,QQ,PolynomialRing,gcd,prod,matrix
out=Path(__file__).resolve().parent;start=time.monotonic();checks={};evidence=[]
raw=load(str(out/'general_curve_zero_gradient_directions.sobj'))
full=load(str(out/'general_curve_dimensionless_full.sobj'))
el=load(str(out/'general_curve_zero_gradient_elimination.sobj'))
ct=load(str(out/'general_curve_zero_gradient_contact.sobj'))
amps=load(str(out/'general_curve_zero_gradient_straight_amplitudes.sobj'))
multi=load(str(out/'general_curve_multiple_zero_gradient.sobj'))
norm=load(str(out/'general_curve_homogeneous_normalized.sobj'))
R=raw['ring'];K=R.base_ring();C=K.ring();H,Z,S,N=K.gens();ts=R.gens()
guards=[C(v) for v in [H,Z,H-Z,Z-1,H+1,H-Z+1,S-1,H*S-Z+1,(H+1)*S-Z]]
def isunit(v,allowed=guards):
    if not v:return False
    for p in [v.numerator(),v.denominator()]:
        for f,e in C(p).factor():
            f=C(f)
            if f.is_constant():continue
            if not any(all((K(f)/g).derivative(x)==0 for x in C.gens()) for g in allowed):return False
    return True
def primitive(v,extra=()):
    p=C(v.numerator())
    if not p:return p
    for g in guards+list(extra):
        while p and not p.quo_rem(g)[1]:p=p.quo_rem(g)[0]
    return C(p/p.leading_coefficient())
# Source scalings are checked in the full original coefficient field.
OC=full['coefficient_ring'];OK=OC.fraction_field();a,r,h,z,s,n=OK.gens()
og=[OC(v) for v in [a,r,1+r,1+a*r,h,z,h-z,z-1,h+1,h-z+1,s-1,h*s-z+1,(h+1)*s-z]]
for i,(f,saved) in enumerate(zip(full['rows'],raw['source_rows'])):
    top={tuple(m):c for m,c in f.dict().items() if sum(m)==(2 if i<3 else 3)};sp={tuple(m):c for m,c in saved.dict().items()}
    key=next(iter(top));ratio=OK(top[key]/sp[key]);good=True
    for p in [ratio.numerator(),ratio.denominator()]:
        for fac,e in OC(p).factor():
            fac=OC(fac)
            if not fac.is_constant() and not any(all((OK(fac)/g).derivative(v)==0 for v in OC.gens()) for g in og):good=False
    checks['full_highest_source_'+str(i)]=set(top)==set(sp) and all(top[m]==ratio*sp[m] for m in top) and good
    checks['spectral_independence_'+str(i)]=all(not c.derivative(a) and not c.derivative(r) for c in sp.values())
expected=[bits for bits in product((0,1),repeat=3) if bits!=(1,1,1)]
checks['seven_patterns']=[tuple(c['nonzero_m']) for c in raw['cases']]==expected
for case in raw['cases']:
    flags=tuple(case['nonzero_m']);tag=''.join(map(str,flags));point=[]
    for j,flag in enumerate(flags):point.extend([ts[j],R.one()] if flag else [R.one(),R.zero()])
    checks[tag+'_direction_normalization']=point==case['point']
    for i,f in enumerate(raw['source_rows']):
        value=R.zero()
        for mon,c in f.dict().items():
            co=K(c.numerator()(0,0,H,Z,S,N))/K(c.denominator()(0,0,H,Z,S,N))
            value+=co*prod(point[j]**mon[j] for j in range(6))
        checks[tag+'_original_direction_source_'+str(i)]=value==case['sources'][i]
print('all seven original direction sources PASS',round(time.monotonic()-start,2),flush=True)

plane_records=[]
for case in el['cases']:
    flags=tuple(case['nonzero_m']);tag=''.join(map(str,flags));pp=next(v['sources'] for v in raw['cases'] if tuple(v['nonzero_m'])==flags)
    active=[j for j,b in enumerate(flags) if b];point=[K.zero()]*3;pivots=[]
    for j,pc in zip(active,case['pivot_cases']):
        idx=pc['source_index'];aa=K(pp[idx].monomial_coefficient(ts[j]));bb=K(pp[idx].constant_coefficient())
        checks[tag+'_original_linear_pivot_'+str(j)]=pp[idx]==aa*ts[j]+bb and aa==pc['a'] and bb==pc['b']
        pa=primitive(aa);pb=primitive(bb);res=C(pa.resultant(pb,C(S)))
        checks[tag+'_pivot_source_strip_'+str(j)]=pa==pc['pivot'] and pb==pc['constant'] and isunit(aa/K(pa)) and isunit(bb/K(pb))
        checks[tag+'_pivot_projection_'+str(j)]=res==pc['resultant'] and bool(res) and not res.degree(C(N)) and not res.degree(C(S))
        for number,(fac,e) in enumerate(res.factor()):
            fac=C(fac)
            if isunit(K(fac)):continue
            if not fac.degree(C(Z)):
                checks[tag+'_pivot_vertical_'+str(j)+'_'+str(number)]=fac.degree(C(H))>0;continue
            den=fac.derivative(C(Z));num=-fac.derivative(C(H));sval=-K(fac.derivative(C(H)))/den
            vals=[K(p(H,Z,sval,0)) for p in [pa,pb]]
            residues=[C(v.numerator()).quo_rem(fac)[1] for v in vals]
            checks[tag+'_pivot_actual_contact_'+str(j)+'_'+str(number)]=fac.gcd(den).is_constant() and gcd([fac]+residues).is_constant()
        point[j]=-bb/aa;pivots.extend(f for f,e in pa.factor())
    checks[tag+'_original_graph']=point==case['point'] and pivots==case['pivots']
    residues=[primitive(K(f(*point)),pivots) for f in pp[3:]]
    checks[tag+'_all_four_trilinear_residues']=residues==case['residuals'] and all(not v.degree(C(N)) for v in residues)
    for i,(f,p) in enumerate(zip(pp[3:],residues)):
        value=K(f(*point));checks[tag+'_residual_only_graph_units_'+str(i)]=isunit(value/K(p),guards+pivots)
    common=gcd(residues);expect=H*S-Z if flags==(0,1,1) else S if flags==(1,0,1) else K.one()
    ratio=K(common)/expect
    checks[tag+'_common_exact']=common==case['common'] and bool(ratio) and all(ratio.derivative(v)==0 for v in C.gens())
    reduced=[p.quo_rem(common)[0] for p in residues];rr=[C(reduced[i].resultant(reduced[j],C(S))) for i,j in combinations(range(4),2)]
    plane=gcd([v for v in rr if v]);quotients=[v.quo_rem(plane)[0] for v in rr if v]
    checks[tag+'_six_resultants']=rr==case['resultants'] and reduced==case['reduced_residuals']
    checks[tag+'_plane_and_finite_remainder']=plane==case['plane'] and all(not v.degree(C(N)) and not v.degree(C(S)) for v in rr) and gcd(quotients).is_constant()
    plane_records.append((flags,plane))
    print('onezero pivots/projections',flags,'PASS',round(time.monotonic()-start,2),flush=True)

# Every nonvertical plane factor is either a saved actual line contact,
# or an inherited guard. Common L and S branches are saved separately.
for flags,plane in plane_records:
    tag=''.join(map(str,flags));saved=[v for v in ct['cases'] if tuple(v['flags'])==flags]
    for j,(fac,e) in enumerate(plane.factor()):
        fac=C(fac)
        if isunit(K(fac)):continue
        if not fac.degree(C(Z)):
            checks[tag+'_plane_vertical_'+str(j)]=fac.degree(C(H))>0;continue
        checks[tag+'_plane_curve_covered_'+str(j)]=any(not K(fac(H,c['Z'],c['S'],0)) for c in saved)
    if flags in ((0,1,1),(1,0,1)):
        checks[tag+'_common_curve_covered']=any((v['S']==Z/H if flags==(0,1,1) else not v['S']) and v['Z']==Z for v in saved)
for i,case in enumerate(ct['cases']):
    flags=tuple(case['flags']);source=next(v['sources'] for v in raw['cases'] if tuple(v['nonzero_m'])==flags)
    graph=next(v['point'] for v in el['cases'] if tuple(v['nonzero_m'])==flags)
    zv=K(case['Z']);sv=K(case['S'])
    def ev(c):return K(c.numerator()(H,zv,sv,0))/K(c.denominator()(H,zv,sv,0))
    values=[ev(K(f(*graph))) for f in source]
    checks['contact_original_values_'+str(i)]=values==case['sources'] and [ev(c) for c in graph]==case['point']
    if zv!=Z:
        checks['line_actual_contact_'+str(i)]=zv.derivative(H)==sv and not sv.derivative(H)
        checks['line_cannot_persist_'+str(i)]=any(values) and all(not v.derivative(Z) and not v.derivative(S) and not v.derivative(N) for v in values)

# The two surviving straight graphs have an entire amplitude row equal
# to a protected coefficient times one squared amplitude, not just one
# nonzero entry in a general row.
for i,case in enumerate(amps['cases']):
    flags=tuple(case['flags']);zv=K(case['Z']);sv=K(case['S']);point=case['point']
    def ev(c):return K(c.numerator()(H,zv,sv,0))/K(c.denominator()(H,zv,sv,0))
    A=matrix(K,6,3)
    for row,f in enumerate(norm['cubics']):
        odd=norm['odd_pairs'][row]
        for mon,c in f.dict().items():
            deg=[mon[2*j]+mon[2*j+1] for j in range(3)];deg[odd]-=1;j=deg.index(2)
            A[row,j]+=ev(c)*prod(K(point[k])**mon[k] for k in range(6))
    checks['straight_amplitude_source_'+str(i)]=A==case['A']
    row,col=(2,0) if flags==(0,1,1) else (0,1)
    expected=-(H-Z)*Z**2/(H**3*(Z-H-1)) if col==0 else -(Z-1)*Z**2/(Z-H-1)
    checks['straight_single_amplitude_contradiction_'+str(i)]=A[row,col]==expected and isunit(expected) and all(not A[row,j] for j in range(3) if j!=col)
checks['both_straight_graphs']=len(amps['cases'])==2

for case in multi['cases']:
    flags=tuple(case['flags']);tag='multi'+''.join(map(str,flags));pp=next(v['sources'] for v in raw['cases'] if tuple(v['nonzero_m'])==flags)
    if not any(flags):
        prim=[primitive(K(p)) for p in pp];base=prim[0];conditions=prim[1:]
    else:
        j=flags.index(1);constants=[primitive(K(p.constant_coefficient())) for p in pp if p.degree()==0 and p]
        checks[tag+'_unique_constant_base']=len(constants)==1;base=constants[0]
        pairs=[(K(p.monomial_coefficient(ts[j])),K(p.constant_coefficient())) for p in pp if p.degree()==1]
        conditions=[primitive(a*d-b*c) for (a,b),(c,d) in combinations(pairs,2) if a*d-b*c]
    checks[tag+'_original_base_and_conditions']=base==case['base'] and conditions==case['conditions'] and not base.degree(C(N)) and all(not p.degree(C(N)) for p in conditions)
    rr=[C(base.resultant(p,C(S))) for p in conditions]
    plane=gcd([v for v in rr if v])
    checks[tag+'_all_resultants']=rr==case['resultants'] and all(not v.degree(C(N)) and not v.degree(C(S)) for v in rr)
    checks[tag+'_finite_plane_intersection']=plane==case['plane'] and bool(plane) and plane.is_constant()
checks['all_multiple_zero_patterns']=len(multi['cases'])==4
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'onezero_plane_data':plane_records,'scope':'Every full-support direction with at least one m_i=0 is excluded on an actual-open ratio curve, including all graph pivots, all seven zero patterns and both surviving straight graphs via original cubic amplitude rows. Resultant gcd is used only for actual curve components; gcd=1 means a finite common plane locus, not pointwise emptiness.'},str(out/'independent_zero_gradient_reader.sobj'))
report={'status':'EXACT_ALL_ZERO_GRADIENT_DIRECTION_PATTERNS_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_zero_gradient_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
