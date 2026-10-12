"""Original-source/candidate-coordinate/canonical-moment determinant replay."""
from pathlib import Path
from itertools import permutations
import time,json
from sage.all import load,save,prod,matrix,vector,QQ

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
original=load(str(out/'general_curve_nonaxis_four_trilinear_probe_input.sobj'))
basis=load(str(out/'independent_two_trilinear_graphs_reader.sobj'))
div=load(str(out/'general_curve_trilinear_candidate_division.sobj'))
claim=load(str(out/'general_curve_trilinear_canonical_matrix.sobj'))
content=load(str(out/'general_curve_trilinear_canonical_content.sobj'))
R=div['ring'];K=R.base_ring();C=K.ring();H,Z,S=K.gens();U,V,T=R.gens();L=H*S-Z;cap=1+H-Z;E3=L+S
oldR=original['ring'];t1,t2,t3=oldR.gens()
images=[V/L-H*T/L+(S-1)/(2*Z),U/S+T/S-(S-1)*L/(2*H*Z*S),T]
checks['exact_inverse_coordinate_images']=images==list(div['images'])
ph=oldR.hom(images,R)
checks['first_proposed_equation_is_U']=ph(basis['proposed_basis'][0])==U
checks['second_proposed_equation_is_V']=ph(basis['proposed_basis'][1])==V
checks['quadratic_is_Q']=ph(basis['proposed_basis'][2])==div['Q']==(2*H*Z*T+L)*(2*Z*cap*T-(S-1)*E3)
inv=R.hom([basis['proposed_basis'][0],basis['proposed_basis'][1],t3],oldR)
for i,v in enumerate(images):checks[f'coordinate_two_sided_inverse_{i}']=inv(v)==oldR.gen(i)
Jac=matrix(R,[[f.derivative(v) for v in (U,V,T)] for f in images])
checks['coordinate_determinant_unit']=Jac.det()==-1/(S*L)
guard=C(H*Z*(H-Z)*(Z-1)*(H+1)*cap*(S-1)*(1+L)*E3*S*L)
denominators={}
for f in images:
    for c in f.coefficients():denominators[str(c.denominator())]=C(c.denominator())
M=claim['matrix'];A=claim['ring'];w=A.gen();checks['canonical_shape']=M.nrows()==M.ncols()==4
def embed(f):return R({(0,0,int(j)):c for j,c in f.dict().items()})
vector=[U*V,U,V,div['Q']]
replayed=[]
for i,f in enumerate(original['sources']):
    transformed=ph(f);it=div['items'][i]
    checks[f'original_full_transformed_{i}']=transformed==it['source'] and f==it['original']
    checks[f'exact_division_source_identity_{i}']=transformed==U*it['A']+V*it['B']+div['Q']*it['C']
    for name in ('B','C'):
        checks[f'no_missed_UV_dependence_{i}_{name}']=not it[name].degree(U) and not it[name].degree(V)
    checks[f'A_affine_only_V_{i}']=not it['A'].degree(U) and it['A'].degree(V)<=1
    canonical=sum((embed(M[i,j])*vector[j] for j in range(4)),R.zero())
    checks[f'complete_canonical_moment_identity_{i}']=canonical==transformed
    replayed.append(transformed)
    for j in range(4):
        for c in M[i,j].coefficients():denominators[str(c.denominator())]=C(c.denominator())
for i,d in enumerate(denominators.values()):
    checks[f'only_inherited_coordinate_and_matrix_denominator_{i}']=all(not guard%C(fac) for fac,e in d.factor())
# Use the Leibniz formula, independently of the creator's determinant routine.
det=A.zero()
for perm in permutations(range(4)):
    sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
    det+=sign*prod(M[i,perm[i]] for i in range(4))
checks['all_24_determinant_terms']=det==claim['determinant']==content['determinant']
checks['det_generic_degree_two']=det.degree()==2
checks['complete_content_identity']=content['primitive']*content['content']==content['denominator']*det
cn=K(content['content'])/(H**6*(Z-H)**2*(Z-1)**2*(H+1)**3)
dn=K(content['denominator'])/(S**2*L**2)
checks['content_exact_old_unit_shape']=cn.numerator().is_constant() and cn.denominator().is_constant() and bool(cn)
checks['denominator_exact_SL_shape']=dn.numerator().is_constant() and dn.denominator().is_constant() and bool(dn)
checks['primitive_three_coefficients_present']=len(content['primitive'].list())==3
for j,c in enumerate(content['primitive'].list()):
    checks[f'primitive_coefficient_polynomial_{j}']=c.denominator()==1
    checks[f'primitive_coefficient_shape_{j}']=c.numerator().total_degree()==16 and len(c.numerator().dict())==[275,238,187][j]
assert all(checks.values()),[k for k,v in checks.items() if not v]
scope='Exact original four-source moment system M(T)*(UV,U,V,Q(T))=0 and determinant/content identity. No determinant coefficient or determinant value is inverted; every zero/degree-drop/identically-zero specialization remains.'
save({'checks':checks,'transformed_sources':replayed,'determinant':det,'content_scalar':cn,'denominator_scalar':dn,'scope':scope},str(out/'independent_trilinear_canonical_reader.sobj'))
summary={'checks_passed':len(checks),'distinct_denominators':len(denominators),'content_scalar':str(cn),'denominator_scalar':str(dn),'elapsed_seconds':time.monotonic()-start,'scope':scope}
(out/'independent_trilinear_canonical_reader.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
