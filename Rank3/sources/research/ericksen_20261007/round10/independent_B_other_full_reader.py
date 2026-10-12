"""Original six-source moments -> other B-root minor -> QQ contact certificate."""
from pathlib import Path
from itertools import permutations
import json,time
from sage.all import load,save,prod,matrix,PolynomialRing,QQ,ZZ,gcd,lcm
out=Path(__file__).resolve().parent; start=time.monotonic(); checks={}
src=load(str(out/'general_curve_centered_six_sources.sobj'))
d=load(str(out/'general_curve_full_six_moment_matrix.sobj'))
mi=load(str(out/'general_curve_full_six_moment_minors.sobj'))
root=load(str(out/'general_curve_B_other_determinant_root.sobj'))
oth=load(str(out/'general_curve_B_other_root_full_minor.sobj'))
ct=load(str(out/'general_curve_B_other_root_contact.sobj'))
A=d['ring']; K=A.base_ring(); C=K.ring(); H,Z,S=K.gens(); L=H*S-Z
R=d['source_ring']; U,V,W=R.gens(); w=A.gen(); images=[(V-H*W)/L,(U+W)/S,W]
checks['coordinate_images']=images==d['images']
X,Y,T=src['ring'].gens(); back=R.hom([S*Y-T,L*X+H*T,T],src['ring'])
checks['coordinate_inverse_all']=all(back(images[i])==src['ring'].gen(i) for i in range(3))
phi=src['ring'].hom(images,R); M=d['matrix']
for i,f in enumerate(src['sources']):
    lifted=phi(f); rebuilt=sum((R(M[i,j](W))*[U*V,U,V,W][j] for j in range(4)),R.zero())
    checks[f'original_full_source_{i}']=lifted==d['sources'][i]==rebuilt
    checks[f'original_constant_divisible_by_W_{i}']=lifted.subs({U:0,V:0,W:0})==0
old=C(H*Z*(H-Z)*(Z-1)*(H+1)*(1+H-Z)*(S-1)*(L+1)*(L+S)*S*L)
denominators={str(c.denominator()):c.denominator() for x in M.list() for c in x.list()}
for i,p in enumerate(denominators.values()): checks[f'protected_source_denominator_{i}']=all(not old%f for f,e in p.factor())
sub=M.matrix_from_rows((0,1,2,3))
minor=sum(((-1)**sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))*prod(sub[i,p[i]] for i in range(4)) for p in permutations(range(4))),A.zero())
checks['independent_24term_original_minor']=minor==mi['items'][0]['determinant']==oth['items'][0]['full_minor']
checks['minor_rows']=tuple(mi['items'][0]['rows'])==tuple(oth['items'][0]['rows'])==(0,1,2,3)
print('six sources / independent determinant PASS',round(time.monotonic()-start,2),flush=True)
F=root['field']; PC=F.ring(); h,z=F.gens(); SP=root['S_ring']; E=root['extension']; ss=E.gen(); B=root['B']
def norm(x):
    c=SP(x.lift()); e=SP((x*ss).lift())
    return matrix(F,[[c[0],e[0]],[c[1],e[1]]]).det()
def ev(p):
    coeff=[F.zero()]*(1+max((int(e[2]) for e in p.dict()),default=0))
    for e,c in p.dict().items(): coeff[int(e[2])]+=F(c)*h**int(e[0])*z**int(e[1])
    value=E.zero()
    for c in reversed(coeff): value=value*ss+c
    return value
inversions={}
def conv(c):
    den=ev(c.denominator()); dn=norm(den); assert dn
    inversions[str(c.denominator())]=dn
    return ev(c.numerator())*den.inverse_of_unit()
ww=root['other_root']-root['simple_root']; checks['centered_other_root']=ww==oth['centered_other_root']
value=E.zero(); coeffs=[]
for j,c in reversed(list(enumerate(minor.list()))):
    cc=conv(c); coeffs.append((j,cc)); value=value*ww+cc
checks['converted_original_coefficients']=coeffs==oth['items'][0]['converted_coefficients']
checks['independent_Horner_value']=value==oth['items'][0]['value']
n=norm(value)
checks['independent_2x2_norm']=n==oth['items'][0]['norm'] and bool(n)
storednorm=oth['items'][0]['norm']
checks['norm_source_contact']=storednorm.numerator()==ct['norm_numerator'] and n.numerator()*storednorm.denominator()==ct['norm_numerator']*n.denominator()
print('B algebra Horner / independent norm PASS',round(time.monotonic()-start,2),flush=True)
f0=ct['norm_numerator']; factor=PC.one()
for i,it in enumerate(ct['stripped_factors']):
    factor*=it['factor']**int(it['exponent'])
    checks[f'stripped_oldguard_{i}']=it['factor'] in [PC(z-1),PC(z-h)]
checks['strip_exact_identity']=f0==factor*ct['after_stripping']
checks['strip_exact_inventory']=[(it['name'],int(it['exponent'])) for it in ct['stripped_factors']]==[('Zminus1',1),('ZminusH',2)]
sq=gcd(ct['after_stripping'],gcd(ct['after_stripping'].derivative(PC.gen(0)),ct['after_stripping'].derivative(PC.gen(1))))
checks['square_gcd_source']=sq==ct['square_gcd'] and sq.is_constant()
f=PC(ct['integer_curve']); fh=f.derivative(PC.gen(0)); fz=f.derivative(PC.gen(1))
ratio=ct['after_stripping']/f
checks['integer_curve_same_zero_set']=ratio.numerator().is_constant() and ratio.denominator().is_constant() and ratio!=0
checks['integral_primitive_coefficients']=all(c.denominator()==1 for c in f.coefficients()) and gcd([ZZ(c) for c in f.coefficients()])==1
HC=PolynomialRing(QQ,'h0')
content=gcd([HC({int(m[0]):c for m,c in f.dict().items() if int(m[1])==j}) for j in range(int(f.degree(PC.gen(1)))+1)])
checks['exact_constant_Z_content']=content.degree()==0 and ct['H_content'].degree()==0
checks['original_curve_derivatives']=fh==ct['F_H'] and fz==ct['F_Z']
checks['full_integer_curve_degree']=f.degree(PC.gen(0))==58 and f.degree(PC.gen(1))==66 and len(f.dict())==2814
cert=ct['certificate']; FF=cert['ring']; k=FF.base_ring(); zz=FF.gen(); hv=int(cert['H_value']); prime=int(cert['prime'])
checks['prime_field']=ZZ(prime).is_prime() and k.characteristic()==prime
def red(p): return FF(sum((k(c)*k(hv)**int(m[0])*zz**int(m[1]) for m,c in p.dict().items()),FF.zero()))
fb,fhb,fzb=red(f),red(fh),red(fz)
checks['full_Z_degree_preserved']=fb.degree()==f.degree(PC.gen(1))==66
checks['integer_curve_reduction']=fb==cert['F_image']
checks['integer_derivative_reductions']=fhb==cert['F_H_image'] and fzb==cert['F_Z_image']
# Reduce exact integer contact coefficient expressions before multiplying;
# this equals reduction of A*F_H^2-D*F_H*F_Z+C*F_Z^2 term by term.
contact=k(hv)*zz*(k(hv)**2-1)*fhb**2-2*(zz**2-k(hv)**2)*fhb*fzb-k(hv)*zz*(zz**2-1)*fzb**2
checks['integer_contact_reduction']=contact==cert['contact_image']
p0,q0=cert['derivative_bezout']; p1,q1=cert['contact_bezout']
checks['derivative_Bezout_exact']=p0*fb+q0*fzb==1
checks['contact_Bezout_exact']=p1*fb+q1*contact==1
inventory={}
for den,dn in inversions.items():
    inventory[den]=[(side,str(p) if p.total_degree()<=2 else f'degree{p.total_degree()}_{len(p.dict())}terms',int(e)) for side,pol in [('numerator',dn.numerator()),('denominator',dn.denominator())] for p,e in pol.factor()]
assert all(checks.values()),[k for k,v in checks.items() if not v]
scope='Exact full-six source and other-root minor replay; degree/content plus integer contact good reduction proves characteristic-zero coprimality. Chart denominator/norm walls are retained and require the separately audited B-normalization, simple-root, old straight-wall, and Mobius exclusions. This is not a whole rank classification.'
summary={'checks_passed':len(checks),'checks':checks,'denominator_norm_inventory':inventory,'elapsed_seconds':time.monotonic()-start,'scope':scope}
save(summary,str(out/'independent_B_other_full_reader.sobj'))
(out/'independent_B_other_full_reader.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
