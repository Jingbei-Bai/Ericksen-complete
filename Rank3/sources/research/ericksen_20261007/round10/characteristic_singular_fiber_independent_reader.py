"""Independent source reader for the singular affine-pair dimension lemma."""
from pathlib import Path
import json,time
from sage.all import load,save,matrix
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
d=load(str(base/'round9/rank1_ratio_characteristic_linear_pair.sobj'))
co=load(str(base/'round9/rank1_ratio_characteristic_coefficients.sobj'))['old_coefficients']
R=d['ring'];u,v,x,s,y=R.gens();K=R.base_ring()
lam,A,B,C,D,E,F,G,J,N=(co[k] for k in ['lam','A','B','C','D','E','F','G','Hcoef','C0'])
alpha=J+lam*D;checks={}
def check(name,condition):checks[name]=bool(condition);assert checks[name],name
f,g=d['linear_source_rows']
M=matrix(R,[[p.derivative(t) for t in (x,y)] for p in (f,g)])
expected=matrix(R,[[u*v+D*s,lam*v*v+B*s*s],[u*u+F*s*s,lam*u*v+J*s]])
check('affine_matrix_from_source',M==expected)
check('affine_no_quadratic_xy',all(not p.derivative(t1).derivative(t2) for p in (f,g) for t1 in (x,y) for t2 in (x,y)))
delta=-B*s*u*u-lam*F*s*v*v+alpha*u*v+s*D*J-B*F*s**3
check('determinant_from_source',M.det()==s*delta)
conic=matrix(R,[[-B*s,alpha/2,0],[alpha/2,-lam*F*s,0],[0,0,s*D*J-B*F*s**3]])
check('homogeneous_conic_determinant',conic.det()==(B*lam*F*s*s-alpha*alpha/4)*s*(D*J-B*F*s*s))
num=d['cramer_numerators'][0]
# Divide the known monomial directly: generic quo_rem is not implemented
# for this Sage multivariate polynomial ring over a function field.
check('cramer_s_divisibility',all(mon[3]>=1 for mon in num.dict()))
X=R({tuple(e-(1 if k==3 else 0) for k,e in enumerate(mon)):c
     for mon,c in num.dict().items()})
check('cramer_s_division_replay',s*X==num)
def uv_coefficient(poly,i,j):
    val=R.zero()
    for mon,c in poly.dict().items():
        check_key=(mon[2]==0 and mon[4]==0)
        assert check_key
        if mon[0]==i and mon[1]==j:val+=c*s**mon[3]
    return val
check('X_u3',uv_coefficient(X,3,0)==0)
check('X_v3',uv_coefficient(X,0,3)==lam*G)
check('X_uv2',uv_coefficient(X,1,2)==lam*(E-A)*s)
check('X_u2v',uv_coefficient(X,2,1)==-lam*C)
check('X_no_uv_degree2',all(not c for mon,c in X.dict().items() if mon[0]+mon[1]==2))
check('candidate_quotient_uv2_obstruction',uv_coefficient(F*s*X+delta*G*v,1,2)==lam*F*(E-A)*s*s+alpha*G)
H=d['remaining_mixed_source'];Jxy=co['Jxy']
check('H_no_x2',not H.derivative(x,2))
check('H_no_y2',not H.derivative(y,2))
check('H_xy_only',H.derivative(x).derivative(y)==Jxy*s)
save({'checks':checks,'source_matrix':M,'delta':delta,'X':X,'H':H,'conic_matrix':conic},str(out/'characteristic_singular_fiber_independent_reader.sobj'))
report={'status':'PASS_INDEPENDENT_SINGULAR_SOURCE_REPLAY','checks':checks,'check_count':len(checks),
 'elapsed_seconds':time.monotonic()-start,
 'scope':'Source identities for the conditional fixed-coefficient lemma. Uniform finite-bad-h factors are reviewed in the independent proof review and the original all-parameter audit.'}
(out/'characteristic_singular_fiber_independent_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'check_count':len(checks),'elapsed_seconds':report['elapsed_seconds']},indent=2),flush=True)
