"""Original four-source input and candidate graphs; no uniform forward lift."""
from pathlib import Path
import json,time
from sage.all import load,save,prod,sage_eval

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
input=load(str(out/'general_curve_nonaxis_four_trilinear_probe_input.sobj'))
unit=load(str(out/'general_curve_unit_gradient_directions.sobj'))
full=load(str(out/'general_curve_dimensionless_full.sobj'))
second=load(str(out/'general_curve_second_trilinear_graph.sobj'))
R=input['ring'];K=R.base_ring();C=K.ring();H,Z,S=K.gens();t1,t2,t3=R.gens();L=H*S-Z;E2=1+L;E3=L+S;cap=1+H-Z
oldK=full['ring'].base_ring();a,r,h,z,s,n=oldK.gens()
def coeff(c):
    assert not c.derivative(a) and not c.derivative(r) and not c.derivative(n)
    return K(c.numerator()(0,0,H,Z,S,0))/K(c.denominator()(0,0,H,Z,S,0))
def conv(f):return R({tuple(m):coeff(c) for m,c in f.dict().items()})
images=[t1,R.one(),t2,R.one(),t3,R.one()]
for i,f in enumerate(full['rows'][3:7]):
    scale=unit['old_cubic_scales'][i]
    rebuilt=sum((coeff(c/scale)*prod(images[j]**int(m[j]) for j in range(6)) for m,c in f.dict().items() if sum(m)==3),R.zero())
    checks[f'original_full_trilinear_source_{i}']=rebuilt==input['sources'][i]==conv(unit['old_cubic_rows'][i])
    guardold=oldK.ring()(a*(a-1)*r*(1+r)*(1+a*r)*h*z*(h-z)*(z-1)*(h+1)*(1+h-z)*(s-1)*(1+h*s-z)*((1+h)*s-z))
    checks[f'trilinear_source_unit_{i}']=all(not guardold%oldK.ring()(fac) for p in (scale.numerator(),scale.denominator()) for fac,e in p.factor())

log=(out/'general_curve_nonaxis_four_trilinear_probe.log').read_text().splitlines();at=log.index('BASIS');symbols=dict(zip(('H','Z','S','t1','t2','t3'),(H,Z,S,t1,t2,t3)))
claimed=[R(sage_eval(v.rstrip(','),locals=symbols)) for v in log[at+1:] if v.strip()]
f1=S*t2-t3+(S-1)*L/(2*H*Z)
f2=L*t1+H*t3-(S-1)*L/(2*Z)
fq=(2*H*Z*t3+L)*(2*Z*cap*t3-(S-1)*E3)
checks['basis_first_literal']=claimed[0]==2*H*Z*f1
checks['basis_second_literal']=claimed[1]==2*Z*f2
checks['basis_third_literal_factorization']=claimed[2]==fq
checks['quadratic_leading_protected']=fq.monomial_coefficient(t3**2)==4*H*Z**2*cap
g0=[S/(2*Z),-L/(2*H*Z),-L/(2*H*Z)]
g1=[-(S-1)*E2/(2*cap*L),(S-1)*E2/(2*H*cap*S),(S-1)*E3/(2*Z*cap)]
checks['second_graph_exact_saved']=g1==list(second['graph'])
guard=C(H*Z*(H-Z)*(Z-1)*(H+1)*cap*(S-1)*E2*E3*S*L)
for label,graph in [('simple',g0),('second',g1)]:
    for j,f in enumerate([f1,f2,fq]):checks[f'{label}_basis_root_{j}']=not f(*graph)
    for j,f in enumerate(input['sources']):checks[f'{label}_original_four_{j}']=not f(*graph)
    for j,g in enumerate(graph):checks[f'{label}_only_protected_denominator_{j}']=all(not guard%C(fac) for fac,e in g.denominator().factor())
P=L**2+H*S**2-Z
checks['collision_polynomial_identity']=P==cap*L+H*(S-1)*E3
checks['first_direction_graph_difference']=g1[0]-g0[0]==-P/(2*Z*cap*L)
checks['second_direction_graph_difference']=g1[1]-g0[1]==P/(2*H*Z*cap*S)
checks['third_direction_graph_difference']=g1[2]-g0[2]==P/(2*H*Z*cap)
mv=[conv(f)(*g1) for f in unit['mixed_differences']]
checks['original_mixed_difference_0']=mv[0]==-(S-1)*P/(H*Z*cap)==second['mixed_differences'][0]
checks['original_mixed_difference_1']=mv[1]==E2*P/(H*Z*cap)==second['mixed_differences'][1]
nv=conv(unit['curvature_value'])(*g1)
checks['second_curvature_from_original_mixed']=nv==(1-S)*E2/(H*cap)==second['N_reconstruction']
checks['curvature_graph_difference']=nv-S*L/(H*Z)==-P/(H*Z*cap)
assert all(checks.values()),[k for k,v in checks.items() if not v]
scope='Original input and two candidate graphs replayed exactly. The second is compatible with mixed equations only at collision and then equals the first. UNIFORM containment of the three proposed basis equations in the original ideal is NOT certified here; all such specialization walls remain.'
save({'checks':checks,'simple_graph':g0,'second_graph':g1,'collision':P,'proposed_basis':[f1,f2,fq],'scope':scope},str(out/'independent_two_trilinear_graphs_reader.sobj'))
summary={'checks_passed':len(checks),'elapsed_seconds':time.monotonic()-start,'scope':scope}
(out/'independent_two_trilinear_graphs_reader.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
