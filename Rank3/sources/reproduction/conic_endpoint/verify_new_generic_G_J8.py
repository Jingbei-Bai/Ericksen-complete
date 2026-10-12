"""Reduce the newly built WHOLE relative G polynomial modulo J8.

No original affine source is used in this reduction: it binds the exact
regular polynomial created by rebuild_generic_G_p21_Q6.py to the dual
conic reference used in the all-orders local-ideal proof.
"""
from sage.all import *
from pathlib import Path
import argparse,json,hashlib,time
from generic_G_decomposition import verify_and_export
ap=argparse.ArgumentParser();ap.add_argument('--source-directory',type=Path,required=True)
ap.add_argument('--output-directory',type=Path,default=Path(__file__).absolute().parent)
args=ap.parse_args();p=args.source_directory;out=args.output_directory;st=time.monotonic()
fresh=load(str(out/'generic_G_relative_regular.sobj'));P=fresh['ring'];PK=P.base_ring();G=fresh['source']
monomial_decomposition=verify_and_export(P,G,out)
export=out/'generic_G_relative_regular_sparse.json'
if not export.exists():
    export.write_text(json.dumps({'coefficient_field':'QQ(k,r,eta,zeta)',
        'variables':list(P.variable_names()),
        'definition':'kappa^13*delta^6*G(p21,Q6) after the exact incidence substitution',
        'terms':[{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in G.dict().items()]},indent=2)+'\n')
bd=load(str(p/'independent_critical_single_RR3_normal.sobj'));BP=bd['ring'];BK=BP.base_ring()
K=FractionField(PolynomialRing(QQ,names=('k','r','eta','zeta','w')));k,r,H,Z,w=K.gens();m=1+k;q=1+r
hom=PK.ring().hom([k,r,H,Z],K);bh=BK.ring().hom([k/r,r,H,Z],K)
c=2*m*(k+w*w)/k;V=PolynomialRing(K,'x');E=V.quotient(V.gen()**2-c,'xb');xb=E.gen()
def basec(a):
    a=BK(a);return E(bh(a.numerator())/bh(a.denominator()))
def scalar(a):
    a=PK(a);return E(hom(a.numerator())/hom(a.denominator()))
def simple(f):
    ans=E.zero()
    for ex,a in f.dict().items():
        assert not any(ex[3:]);ans+=basec(a)*xb**ex[0]*(-xb)**ex[1]*w**ex[2]
    return ans
C=bd['graph_boundary'][27];ox,oz=BP.gens()[:2]
oa,orr,oH,oZ=BK.gens();ow=BP.gen(2);ok=oa*orr;om=1+ok
explicit_C=2*ok/om**2*(ox**3*oz+ox*oz**3)+(10*ox**2-4*ox*oz-6*oz**2)*ow**2/om+ok*(-6*ox**2-4*ox*oz+10*oz**2)/om
assert C==explicit_C,'stored C does not equal the freshly reconstructed displayed polynomial'
Cx=simple(C.derivative(ox));Cz=simple(C.derivative(oz))
xs=E((3*w*w-k)/(2*(k+w*w)))
assert (Cx-Cz)*xs+Cz==0
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def mul(a,b):return (a[0]*b[0],a[0]*b[1]+a[1]*b[0])
def scale(a,b):return (a[0]*b,a[1]*b)
one=(E.one(),E.zero());zero=(E.zero(),E.zero());ss=(E.zero(),E.one())
xx=(xb,xs);zz=(-xb,1-xs);ww=(E(w),E.zero())
f0=add(scale(mul(xx,ww),-2),scale(mul(ss,ww),(4*H*Z-2*H-4*Z**2+2*Z)/(H*Z)))
f1=scale(ww,-2*Z*(H-1)/(H*(Z-1)))
f2=add(scale(xx,2),scale(ss,(2*H*Z-2*H-4*Z**2+2*Z)/Z))
dd=scale(mul(ss,ww),4*r*(H+1)**3/(q*H**2))
images=(xx,zz,ww,f0,f1,f2,dd,zero,zero)
cache={}
def power(j,n):
    if (j,n) not in cache:
        a=one
        for _ in range(n):a=mul(a,images[j])
        cache[j,n]=a
    return cache[j,n]
ans=zero;kept=0
for ex,a in G.dict().items():
    if ex[7] or ex[8] or ex[6]>=2:continue
    term=one
    for j,n in enumerate(ex):
        if n:term=mul(term,power(j,int(n)))
    ans=add(ans,scale(term,scalar(a)));kept+=1
assert ans==zero,('new complete G is not in J8',ans)
# Verify the defining reference C itself, not an assumed conic equation.
Cdual=zero
for ex,a in C.dict().items():
    assert not any(ex[3:]);term=one
    for j,n in enumerate(ex[:3]):
        if n:term=mul(term,power(j,int(n)))
    Cdual=add(Cdual,scale(term,basec(a)))
assert Cdual==zero
# The dual lift extends across w=0: every denominator coefficient is
# nonzero after w=0,eta=1,zeta=-1, with k,r indeterminate.
K0=FractionField(PolynomialRing(QQ,names=('kk','rr')));kk,rr=K0.gens()
cyc=K.ring().hom([kk,rr,1,-1,0],K0)
denoms=[]
for coord in [xs]+[a for pair in images for a in pair]:
    for a in E(coord).list():
        val=cyc(K(a).denominator());assert val!=0
        denoms.append(str(val))
report={'status':'PASS_NEW_GENERIC_G_EXACT_J8_MEMBERSHIP','scope':__doc__,
        'generic_source_sha256':hashlib.sha256((out/'generic_G_relative_regular.sobj').read_bytes()).hexdigest(),
        'full_source_monomials':len(G.dict()),'surviving_reference_monomials':kept,
        'conic_equation':'xb^2=2*(1+k)*(k+w^2)/k','dual_lift_x_s':str(xs),
        'C_equals_displayed_freshly_verified_C0':True,
        'dual_lift_equals_(3w2-k)/(2(k+w2))':True,
        'complete_G_membership_in_t_kappa':monomial_decomposition,
        'C_constant_and_dual_coefficients_zero':True,'G_constant_and_dual_coefficients_zero':True,
        'dual_lift_regular_at_w_zero_for_all_allowed_materials':True,
        'seconds':time.monotonic()-st}
(out/'new_generic_G_J8.json').write_text(json.dumps(report,indent=2)+'\n')
save({'ring':E,'G_remainder':ans,'C_remainder':Cdual,'reference_images':images,'report':report},
     str(out/'new_generic_G_J8.sobj'),compress=False)
print(json.dumps(report,indent=2),flush=True)
