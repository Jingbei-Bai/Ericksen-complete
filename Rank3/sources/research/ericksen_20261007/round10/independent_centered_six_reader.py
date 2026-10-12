"""Full old-source bridge and all twenty centered Jacobian minors."""
from pathlib import Path
from itertools import combinations
import json,time
from sage.all import load,save,prod,matrix,PolynomialRing,sage_eval,gcd

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
d=load(str(out/'general_curve_centered_six_sources.sobj'))
unit=load(str(out/'general_curve_unit_gradient_directions.sobj'))
full=load(str(out/'general_curve_dimensionless_full.sobj'))
R=d['ring'];K=R.base_ring();C=K.ring();H,Z,S=K.gens();X,Y,T=R.gens();L=H*S-Z;P=L**2+H*S**2-Z
B=H*Z*(H**2-1)*S**2+2*(Z**2-H**2)*S-H*Z*(Z**2-1)
oldK=full['ring'].base_ring();a,r,h,z,s,n=oldK.gens()
OR=PolynomialRing(oldK,names=('t1','t2','t3'));ts=OR.gens();images=[ts[0],OR.one(),ts[1],OR.one(),ts[2],OR.one()]
def top_eval(f,deg):return sum((c*prod(images[j]**int(m[j]) for j in range(6)) for m,c in f.dict().items() if sum(m)==deg),OR.zero())
mix=[top_eval(f,2)*r/(1+r) for f in full['rows'][:3]]
source=[mix[1]-mix[0],mix[2]-mix[0]]+[top_eval(f,3)/unit['old_cubic_scales'][i] for i,f in enumerate(full['rows'][3:7])]
center=[X+S/(2*Z),Y-L/(2*H*Z),T-L/(2*H*Z)]
checks['complete_center_coordinates']=center==list(d['simple_center'])
checks['all_six_sources_present']=len(source)==len(d['sources'])==6
def cv(c):
    assert not c.derivative(a) and not c.derivative(r) and not c.derivative(n)
    return K(c.numerator()(0,0,H,Z,S,0))/K(c.denominator()(0,0,H,Z,S,0))
for i,f in enumerate(source):
    rebuilt=sum((cv(c)*prod(center[j]**int(m[j]) for j in range(3)) for m,c in f.dict().items()),R.zero())
    checks[f'original_full_source_centered_{i}']=rebuilt==d['sources'][i]
    checks[f'center_is_exact_zero_{i}']=rebuilt(0,0,0)==0
J=matrix(K,[[f.derivative(v)(0,0,0) for v in R.gens()] for f in d['sources']])
checks['all_eighteen_Jacobian_entries']=J==d['linear_matrix']
expected=list(combinations(range(6),3));checks['all_twenty_row_triples']=len(d['records'])==20 and [tuple(v['rows']) for v in d['records']]==expected
minors={};triquot={};symbols={'H':H,'Z':Z,'S':S}
for item in d['records']:
    rows=tuple(item['rows']);value=J.matrix_from_rows(rows).det();minors[rows]=value
    checks['minor_zero_flag_'+str(rows)]=(not bool(value))==item['zero']
    parsed=K.zero() if item['factorization'] is None else K(sage_eval(item['factorization'],locals=symbols))
    checks['complete_minor_factorization_'+str(rows)]=value==parsed
    if all(i>=2 for i in rows):
        quo=value/(P*B);triquot[rows]=quo
        checks['old_four_PB_divides_numerator_'+str(rows)]=not C(value.numerator())%C(P*B)
guard=C(H*Z*(H-Z)*(Z-1)*(H+1)*(1+H-Z)*(S-1)*(1+H*S-Z)*((1+H)*S-Z)*S*L)
guards=[C(fac) for fac,e in guard.factor()]
def strip(p):
    p=C(p)
    for fac in guards:
        while p and not p%fac:p=p//fac
    return p
residuals=[strip(q.numerator()) for q in triquot.values()]
common=gcd(residuals)
checks['old_four_have_exactly_PB_common_nonunit_factor']=common.is_constant() and bool(common)
checks['B_agrees_with_contact_source']=B==load(str(out/'general_curve_B_wall_contact.sobj'))['B']
assert all(checks.values()),[k for k,v in checks.items() if not v]
scope='Exact full-source centering and all 20 origin-Jacobian minors. Old-four trilinear minors have exactly common P*B after inherited units; this does not give global direction coverage or exclude other components.'
save({'checks':checks,'jacobian':J,'minors':minors,'P':P,'B':B,'old_four_quotients':triquot,'old_four_residual_gcd':common,'scope':scope},str(out/'independent_centered_six_reader.sobj'))
summary={'checks_passed':len(checks),'old_four_P':str(P),'old_four_B':str(B),'elapsed_seconds':time.monotonic()-start,'scope':scope}
(out/'independent_centered_six_reader.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
