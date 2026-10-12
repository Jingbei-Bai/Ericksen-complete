"""Quadratic-algebra norms and exact characteristic-zero good reduction."""
from pathlib import Path
import json,time
from sage.all import load,save,prod,matrix,PolynomialRing,QQ,ZZ,gcd,lcm

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
d=load(str(out/'general_curve_B_other_determinant_root.sobj'))
src=d['source'];orig=load(str(out/'general_curve_trilinear_canonical_matrix.sobj'))
certs=load(str(out/'general_curve_B_norm_contact_certificates.sobj'))
F=d['field'];C=F.ring();H,Z=F.gens();sp=d['S_ring'];s=sp.gen();E=d['extension'];ss=E.gen();W=src['ring'];w=W.gen()
B=H*Z*(H**2-1)*s**2+2*(Z**2-H**2)*s-H*Z*(Z**2-1)
checks['original_B_polynomial']=B==src['B'] and B.monic()==d['B']
checks['quadratic_algebra_relation']=E(B)==0 and d['B'].degree()==2
def norm(v):
    p=sp(v.lift());p2=sp((v*ss).lift())
    return matrix(F,[[p[0],p2[0]],[p[1],p2[1]]]).det()
inversions={}
def peval(p):return sum((E(c)*E(H)**int(m[0])*E(Z)**int(m[1])*ss**int(m[2]) for m,c in p.dict().items()),E.zero())
def coeff(c):
    den=peval(c.denominator());dn=norm(den);assert dn
    inversions[str(c.denominator())]=dn
    return peval(c.numerator())*den.inverse_of_unit()
replayed=W([coeff(c) for c in orig['determinant'].list()])
checks['original_canonical_det_reduced']=replayed==src['determinant']
checks['simple_root_formula']=d['simple_root']==-(H*ss-Z)/(2*H*Z)
checks['simple_is_original_det_root']=replayed(d['simple_root'])==0
checks['leading_coefficient_source']=d['leading']==replayed[2]
checks['leading_norm_independent_matrix']=norm(d['leading'])==d['leading_norm']
checks['difference_norm_independent_matrix']=norm(d['other_root']-d['simple_root'])==d['difference_norm']
checks['other_root_exact_formula']=d['leading']*(d['other_root']+d['simple_root'])==-replayed[1]
checks['other_root_original_det_zero']=replayed(d['other_root'])==0
checks['all_coefficient_root_factorization']=replayed==d['leading']*(w-d['simple_root'])*(w-d['other_root'])
print('quadratic norm/source replay PASS',round(time.monotonic()-start,2),flush=True)

norms={'leading_norm':d['leading_norm'],'root_difference_norm':d['difference_norm']}
factors={key:list(value.numerator().factor()) for key,value in norms.items()}
expected=[]
for key,ff in factors.items():
    for pol,exponent in ff:
        if pol.degree(C.gen(1)) and pol.total_degree()>1:
            primitive=C(pol/pol.leading_coefficient())
            if not any(primitive==f for k,f in expected):expected.append((key,primitive))
checks['complete_nontrivial_norm_curve_list']=len(expected)==len(certs['items'])==3
records=[]
for idx,item in enumerate(certs['items']):
    f=C(item['curve']);key=item['source'];expectedf=next(pp for kk,pp in expected if kk==key and f/f.leading_coefficient()==pp)
    checks[f'norm_factor_source_{idx}']=not C(norms[key].numerator())%f
    checks[f'integer_curve_coefficients_{idx}']=all(c.denominator()==1 for c in f.coefficients())
    checks[f'primitive_integer_content_{idx}']=gcd([ZZ(c) for c in f.coefficients()])==1
    HC=PolynomialRing(QQ,'hh0');hh0=HC.gen()
    zcontent=gcd([HC({int(m[0]):c for m,c in f.dict().items() if int(m[1])==j}) for j in range(int(f.degree(C.gen(1)))+1)])
    checks[f'exact_constant_Z_content_{idx}']=zcontent.degree()==0 and item['Z_content'].degree()==0 and QQ(zcontent[0])==QQ(item['Z_content'][0])
    fh=f.derivative(C.gen(0));fz=f.derivative(C.gen(1))
    checks[f'original_derivatives_{idx}']=fh==item['F_H'] and fz==item['F_Z']
    contact=C(H*Z*(H**2-1))*fh**2-2*C(Z**2-H**2)*fh*fz-C(H*Z*(Z**2-1))*fz**2
    cert=item['certificate']
    if cert is None:
        checks['small_curve_kept_UNRESOLVED']=f==C(H*Z-H+Z)
        checks['small_curve_is_B_integral']=not contact%f
        checks['small_curve_explicit_B_root']=B.subs({s:1/(H+1)**2}).subs({Z:H/(H+1)})==0
        records.append({'index':idx,'status':'INTEGRAL_BASE_CURVE_RETAINED','curve':str(f)})
        continue
    A=cert['ring'];k=A.base_ring();zz=A.gen();prime=int(cert['prime']);hv=int(cert['H_value'])
    checks[f'prime_field_{idx}']=int(k.characteristic())==prime and ZZ(prime).is_prime()
    def ev(pol):return sum((k(c)*k(hv)**int(m[0])*zz**int(m[1]) for m,c in pol.dict().items()),A.zero())
    fb=ev(f);fhb=ev(fh);fzb=ev(fz);ct=ev(contact)
    checks[f'exact_integer_curve_reduction_{idx}']=fb==cert['F_image']
    checks[f'exact_derivative_reductions_{idx}']=fhb==cert['F_H_image'] and fzb==cert['F_Z_image']
    checks[f'exact_integer_contact_reduction_{idx}']=ct==cert['contact_image']
    checks[f'preserved_FULL_Z_degree_{idx}']=fb.degree()==f.degree(C.gen(1))==14
    u0,v0=cert['derivative_bezout'];u1,v1=cert['contact_bezout']
    checks[f'derivative_Bezout_full_replay_{idx}']=u0*fb+v0*fzb==1
    checks[f'contact_Bezout_full_replay_{idx}']=u1*fb+v1*ct==1
    records.append({'index':idx,'total_degree':int(f.total_degree()),'degree_Z':int(f.degree(C.gen(1))),'terms':len(f.dict()),'prime':prime,'H_value':hv,'status':'EXACT_QQ_COPRIME_BY_DEGREE_AND_CONTENT'})

# Enumerate, rather than discard, all other norm factors.
inventory={}
for key,val in norms.items():
    inventory[key]=[(side,str(pol) if pol.total_degree()<=2 else f'degree{pol.total_degree()}_{len(pol.dict())}terms',int(exponent)) for side,p in [('numerator',val.numerator()),('denominator',val.denominator())] for pol,exponent in p.factor()]
assert all(checks.values()),[k for k,v in checks.items() if not v]
scope='Exact B-algebra roots and two QQ plane-contact exclusions from integer good reduction. F=HZ-H+Z is a genuine integral base curve retained by this contact certificate; its separate full-direction exclusion is independent_B_mobius_unit_reader. Norm inversion is only valid off its explicitly retained zero/denominator loci.'
save({'checks':checks,'records':records,'norm_factor_inventory':inventory,'inverted_original_denominator_norms':inversions,'scope':scope},str(out/'independent_B_norm_contact_reader.sobj'))
summary={'checks_passed':len(checks),'records':records,'norm_factor_inventory':inventory,'elapsed_seconds':time.monotonic()-start,'scope':scope}
(out/'independent_B_norm_contact_reader.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
