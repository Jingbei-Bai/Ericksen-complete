"""Independent all-spectrum scaling, leading-source and contact audit.

Only uniquely named *_full and *_leading artifacts are used.  Scalar
contact equations are conditional on an actual open scalar-monodromy wall.
"""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing,matrix

out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
old=load(str(out/'general_curve_full_gradient_sources.sobj'))
full=load(str(out/'general_curve_dimensionless_full.sobj'))
lead=load(str(out/'general_curve_dimensionless_leading.sobj'))
tri=load(str(out/'general_curve_mixed_triangle.sobj'))
sc=load(str(out/'general_curve_scalar_dimensionless.sobj'))
contact=load(str(out/'general_curve_scalar_contact_source.sobj'))
R=full['ring'];K=R.base_ring();CR=K.ring();a,r,H,Z,S,N=K.gens();q=1+r
xs=R.gens();args=(a,r,H/r,q*Z/r,q*S,q*r*N)
def cc(c):return K(c.numerator()(*args))/K(c.denominator()(*args))
def convert(p):
    return R({tuple(mm):cc(c)/r**sum(mm[j] for j in (1,3,5)) for mm,c in p.dict().items()})
sources=[old['rows'][ij] for ij in [(0,1),(0,2),(1,2)]]+old['cubic_sources']
for i,p in enumerate(sources):
    checks['full_source_'+str(i)]=convert(p)==full['rows'][i]
    checks['full_primitive_'+str(i)]=full['rows'][i]==full['primitive_rows'][i]*full['scalings'][i]

allowed=[a,r,q,1+a*r,H,Z,Z-H,1-Z,1+H,1+H-Z,1-S,1+H*S-Z,(1+H)*S-Z]
def protected_unit(c):
    c=K(c)
    for poly in [c.numerator(),c.denominator()]:
        for factor,power in poly.factor():
            if factor.is_constant():continue
            if not any(all((K(factor)/guard).derivative(x)==0 for x in K.gens()) for guard in allowed):
                return False
    return bool(c)
for i,scale in enumerate(full['scalings']):checks['full_scaling_protected_'+str(i)]=protected_unit(scale)
def lk(poly):return K(poly(*K.gens()))
for pos,ij in enumerate([(0,1),(0,2),(1,2)]):
    B=lead['blocks'][ij];num,den=lead['block_scales'][ij];scale=lk(num)/lk(den)
    checks['mixed_scaling_protected_'+str(pos)]=protected_unit(scale)
    for i in range(2):
        for j in range(2):
            c=full['rows'][pos].monomial_coefficient(xs[2*ij[0]+i]*xs[2*ij[1]+j])
            checks['mixed_leading_'+str(pos)+str(i)+str(j)]=c==scale*lk(B[i,j])
            checks['mixed_spectrum_free_'+str(pos)+str(i)+str(j)]=lk(B[i,j]).derivative(a)==0 and lk(B[i,j]).derivative(r)==0
for k,row in enumerate(lead['cubic_rows']):
    num,den=lead['cubic_scales'][k];scale=lk(num)/lk(den)
    checks['cubic_scaling_protected_'+str(k)]=protected_unit(scale)
    for col,(i,j,l) in enumerate((i,j,l) for i in range(2) for j in range(2) for l in range(2)):
        c=full['rows'][k+3].monomial_coefficient(xs[i]*xs[2+j]*xs[4+l])
        checks['cubic_leading_'+str(k)+'_'+str(col)]=c==scale*lk(row[col])
        checks['cubic_spectrum_free_'+str(k)+'_'+str(col)]=lk(row[col]).derivative(a)==0 and lk(row[col]).derivative(r)==0

# Guard conversion includes every noncharacteristic slope and K, with no
# assumption that h,z,kappa,nu or gradient components are positive.
oldguards=[H/r,q*Z/r,q*(Z-H)/r,1-Z,1+H,1+H-Z,q*(1-S),q*(1+H*S-Z),q*((1+H)*S-Z)]
newguards=[H,Z,Z-H,1-Z,1+H,1+H-Z,1-S,1+H*S-Z,(1+H)*S-Z]
for i,(o,n) in enumerate(zip(oldguards,newguards)):checks['guard_ratio_'+str(i)]=protected_unit(o/n)

# Re-evaluate the three ORIGINAL rational scalar defects, not only their
# normalized numerators. Coefficients of saved denominator polynomials lie
# in QQ(r), so embed that ring explicitly instead of relying on names.
def from_r(c):return K(c.numerator()(r))/K(c.denominator()(r))
def from_scalar_poly(p):
    return sum((from_r(c)*H**mm[0]*Z**mm[1]*S**mm[2]*N**mm[3] for mm,c in p.dict().items()),K.zero())
for i,e in enumerate(tri['scalar_defects']):
    p=K(sc['polynomials'][i](H,Z,S,N));pref=from_r(sc['prefactors_in_r'][i]);den=from_scalar_poly(sc['denominators'][i])
    checks['original_scalar_defect_'+str(i)]=cc(e)==pref*p/den
    checks['scalar_prefactor_protected_'+str(i)]=protected_unit(pref)
    checks['scalar_denominator_protected_'+str(i)]=protected_unit(den)

C=contact['ring'];HH,ZZ,SS,NN,WW=C.gens();embed=sc['ring'].hom([HH,ZZ,SS,NN],C)
for i,p in enumerate(sc['polynomials']):
    row=contact['rows'][i];mul=contact['removed_guard_content'][i]
    checks['removed_guard_identity_'+str(i)]=embed(p)==mul*row
    checks['removed_guard_exact_'+str(i)]=mul==[C.one(),1+HH-ZZ,C.one()][i]
    # Independent monomial-chain evaluation of d/dH.
    derivative=C.zero()
    for mm,c in row.dict().items():
        for j,velocity in enumerate([C.one(),SS,NN,WW]):
            if mm[j]:
                shifted=list(mm);shifted[j]-=1
                derivative+=c*mm[j]*C.monomial(*shifted)*velocity
    checks['contact_derivative_'+str(i)]=derivative==contact['contact_derivatives'][i]
checks['nine_contact_guards']=contact['guards']==[HH,ZZ,ZZ-HH,1-ZZ,1+HH,1+HH-ZZ,1-SS,1+HH*SS-ZZ,(1+HH)*SS-ZZ]

# Algebraic chain rule in original curve variables (h,z,kappa,nu,nu_h).
J=PolynomialRing(QQ,names=('r','h','z','kappa','nu','nu_h'));F=J.fraction_field();rr,hh,zz,kk,nn,nn_h=F.gens();qq=1+rr
def dh(f):return f.derivative(hh)+kk*f.derivative(zz)+nn*f.derivative(kk)+nn_h*f.derivative(nn)
checks['contact_dZ_dH']=dh(rr*zz/qq)/rr==kk/qq
checks['contact_dS_dH']=dh(kk/qq)/rr==nn/(qq*rr)
checks['contact_dN_dH']=dh(nn/(qq*rr))/rr==nn_h/(qq*rr**2)

assert all(checks.values()),[k for k,v in checks.items() if not v]
payload={'checks':checks,'allowed_factors':allowed,'source_files':[
    'general_curve_full_gradient_sources.sobj','general_curve_dimensionless_full.sobj',
    'general_curve_dimensionless_leading.sobj','general_curve_mixed_triangle.sobj',
    'general_curve_scalar_dimensionless.sobj','general_curve_scalar_contact_source.sobj'],
    'scope':'All-spectrum exact source/scaling/contact audit. Highest parts alone are spectrum-free; full lower-degree sources retain spectral parameters. Contact constraints are conditional on scalar-monodromy holding on an actual open curve interval.'}
save(payload,str(out/'independent_dimensionless_audit.sobj'))
report={'status':'EXACT_ALL_PARAMETER_AUDIT_PASS','check_count':len(checks),'checks':checks,
        'elapsed_seconds':time.monotonic()-start,'scope':payload['scope']}
(out/'independent_dimensionless_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
