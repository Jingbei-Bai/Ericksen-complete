"""Independent replay: full simple-T moment fiber on an actual B interval."""
from pathlib import Path
from itertools import combinations,permutations
import json,time
from sage.all import load,save,matrix,prod,gcd,QQ

out=Path(__file__).resolve().parent; start=time.monotonic(); checks={}
can=load(str(out/'general_curve_trilinear_canonical_matrix.sobj'))
mom=load(str(out/'general_curve_B_simpleT_moments.sobj'))
fac=load(str(out/'general_curve_B_simpleT_minor_factors.sobj'))
pr=load(str(out/'general_curve_B_simpleT_projection.sobj'))
ct=load(str(out/'general_curve_B_simpleT_contact_exclusion.sobj'))
K=mom['original_field']; C=K.ring(); H,Z,S=K.gens(); L=H*S-Z
alpha=-L/(2*H*Z)
M=matrix(K,[[can['matrix'][i,j](alpha) for j in range(3)] for i in range(4)])
checks['simple_root_exact']=mom['T_root']==alpha
div=load(str(out/'general_curve_trilinear_candidate_division.sobj'))
checks['original_canonical_Q_zero']=div['Q'].subs({div['ring'].gen(2):alpha})==0
for i in range(4):
    for j in range(3): checks[f'original_matrix_{i}_{j}']=M[i,j]==mom['matrix'][i,j]
old=C(H*Z*(H-Z)*(Z-1)*(H+1)*(1+H-Z)*(S-1)*(L+1)*(L+S)*S*L)
for i,v in enumerate(M.list()):
    checks[f'protected_matrix_denominator_{i}']=all(not old%f for f,e in v.denominator().factor())
R=pr['ring']; F=R.base_ring(); B=pr['B']; h,z=F.gens(); s=R.gen(); PC=F.ring()
def tr(p): return R(sum((F(c)*h**int(m[0])*z**int(m[1])*s**int(m[2]) for m,c in p.dict().items()),R.zero()))
def det3(a):
    return sum(((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))*prod(a[i,p[i]] for i in range(3)) for p in permutations(range(3))),K.zero())
checks['B_original_formula']=B==h*z*(h*h-1)*s*s+2*(z*z-h*h)*s-h*z*(z*z-1)
checks['monic_B_exact']=pr['monic_B']==B.monic()
minor_polys=[]; norms=[]
for k,ids in enumerate(combinations(range(4),3)):
    v=det3(M.matrix_from_rows(ids)); stored=mom['items'][k]['original']; p=stored.numerator(); minor_polys.append(p)
    checks[f'full_original_minor_{k}']=v==stored and v.numerator()*stored.denominator()==p*v.denominator()
    checks[f'original_minor_remainder_{k}']=tr(p)%B==mom['items'][k]['remainder']
    checks[f'exact_common_factor_division_{k}']=p==fac['common_numerator_factor']*fac['quotients'][k]
    item=pr['items'][k]; f=tr(fac['quotients'][k]); rem=item['remainder']; quo=item['B_division_quotient']
    checks[f'projection_source_{k}']=f==item['source_in_S'] and fac['quotients'][k]==item['original_quotient']
    checks[f'polynomial_B_division_{k}']=f==quo*B+rem and rem.degree()<=1
    u,vv=rem[0],rem[1]; b=B.monic()
    # Independent norm as a multiplication matrix in F[s]/(B).
    nv=matrix(F,[[u,-b[0]*vv],[vv,u-b[1]*vv]]).det(); norms.append(PC(nv.numerator()))
    checks[f'independent_norm_matrix_{k}']=nv==item['norm'] and PC(nv.numerator())==item['norm_numerator']
    checks[f'norm_identity_{k}']=rem*(u-b[1]*vv-vv*s)==nv-vv*vv*b
    gb=PC(h*z*(h*h-1))
    checks[f'finite_H_only_norm_denominators_{k}']=all(not gb%ff for ff,e in nv.denominator().factor())
checks['actual_common_minor_gcd']=gcd(minor_polys)==fac['common_numerator_factor']
cf=fac['common_numerator_factor']
checks['common_minor_factor_oldunit']=all(not old%ff for ff,e in cf.factor())
g=PC(ct['norm_common_factor']); curve=PC(ct['rankdrop_curve'])
checks['common_all_norms']=gcd(norms)==pr['common_norm_numerator']==g
checks['first_two_norms_suffice']=gcd(norms[:2])==g
checks['residual_norms_exact']=[n//g for n in norms[:2]]==ct['norm_residuals']
checks['residual_norms_coprime']=gcd([n//g for n in norms[:2]]).is_constant()
checks['only_oldguard_factors_removed_from_projection']=all(not PC(h*z*(h-z)*(z-1)*(h+1))%ff for ff,e in (g//curve).factor())
ch,cz=curve.derivative(PC.gen(0)),curve.derivative(PC.gen(1))
J=PC(B[2])*ch*ch-PC(B[1])*ch*cz+PC(B[0])*cz*cz
checks['rankdrop_curve_degrees']=curve.total_degree()==9 and len(curve.dict())==24
checks['curve_derivatives_exact']=ct['curve_derivatives']==[ch,cz]
checks['full_contact_polynomial_source']=J==ct['contact_projection']
checks['contact_identity_without_RZ_division']=R(J)==B[2]*(R(ch)+s*R(cz))**2-(2*B[2]*s+B[1])*R(cz)*(R(ch)+s*R(cz))+B*R(cz)**2
checks['contact_exact_QQ_coprime']=gcd(curve,J).is_constant()
checks['contact_size']=J.total_degree()==20 and len(J.dict())==104
assert all(checks.values()),[k for k,v in checks.items() if not v]
scope='The entire simple-T moment branch on a persistent real B interval: all U,V and rank-drop cases retained. Plane coprimality excludes persistent rank-drop, so dense full rank forces the simple graph, already excluded on B. No pointwise absence or global rank classification is claimed.'
summary={'checks_passed':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start,'scope':scope}
save(summary,str(out/'independent_B_simpleT_reader.sobj'))
(out/'independent_B_simpleT_reader.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
