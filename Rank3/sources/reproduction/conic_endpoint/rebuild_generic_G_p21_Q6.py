"""Build the COMPLETE relative G(p21,Q6) before any conic quotient.

Parameters k,r,eta,zeta and all nine incidence variables are independent.
Every negative kappa Newton layer is proved polynomial in the WHOLE
incidence ring; no restriction to J8 is used to obtain regularity.
"""
from sage.all import *
from pathlib import Path
import argparse,hashlib,json,time
from generic_G_decomposition import verify_and_export
ap=argparse.ArgumentParser();ap.add_argument('--source-directory',type=Path,required=True)
ap.add_argument('--output-directory',type=Path,default=Path(__file__).absolute().parent)
args=ap.parse_args();p=args.source_directory;out=args.output_directory;out.mkdir(parents=True,exist_ok=True)
st=time.monotonic();gn='independent_partial_D_Newton_gauges_reader.sobj';en='independent_delta_subcritical_even_full_source.sobj'
g=load(str(p/gn));ev=load(str(p/en));old=g['ring'].base_ring();OK=ev['ring'].base_ring()
oa,orr,oH,oZ=OK.gens()
K=FractionField(PolynomialRing(QQ,names=('k','r','eta','zeta')));k,r,H,Z=K.gens();m=1+k;q=1+r
ph=old.ring().hom([k/r,k*q/(r*m),H/r,Z*q/r],K);po=OK.ring().hom([k/r,r,H,Z],K)
R=PolynomialRing(K,names=ev['ring'].variable_names());vs=R.gens();cache={}
def cp(f):
    field=f.parent().base_ring();hh=ph if field==old else po;result={}
    for ex,c in f.dict().items():
        key=(str(field),c)
        if key not in cache:cache[key]=hh(c.numerator())/hh(c.denominator())
        result[tuple(map(int,ex))]=cache[key]
    return R(result)
gi=2;qi=6;it=g['items'][gi]
fv=[cp(f) for f in it['field_numerators']];bv=[cp(f) for f in it['base_numerators']]
OQ=next(v['full_source'] for v in ev['items'] if int(v['Q_index'])==qi);Q=cp(OQ)
G=R.zero()
for j in range(9):
    G+=fv[j]*Q.derivative(vs[j])
    print('derivative summand',j,'terms',len(G.dict()),'seconds',time.monotonic()-st,flush=True)
for j,param in enumerate((oH,oZ)):
    scale=orr if j==0 else orr/(1+orr)
    derivative=OQ.parent()({tuple(map(int,e)):scale*c.derivative(param)
                           for e,c in OQ.dict().items() if c.derivative(param)})
    G+=bv[j]*cp(derivative)
    print('base derivative',j,'terms',len(G.dict()),'seconds',time.monotonic()-st,flush=True)
save({'ring':R,'G':G,'scope':__doc__},str(out/'generic_G_original_polynomial.sobj'),compress=False)
P=PolynomialRing(K,names=('x','z','w','e0','e1','e2','delta','t','kappa'))
x,z,w,e0,e1,e2,delta,t,kap=P.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
images=(x,P.one(),z,w,(t*z-kap*e2)/beta,x*w*t/beta+kap*(x*e1-e0)/alpha,
        alpha*w*t/beta+kap*e1,t,P.one())
hom=R.hom(images,P);layers={}
for ex,c in G.dict().items():
    degree=sum(ex);rdegree=sum(ex[6:]);layers.setdefault((degree,rdegree),{})[tuple(map(int,ex))]=c
ans=P.zero();records=[];layer_polys=[]
for (degree,rdegree),coeff in sorted(layers.items()):
    kp=13-2*degree-rdegree;dp=6-degree;assert dp>=0
    pulled=hom(R(coeff));before=len(pulled.dict())
    min_kappa=min((int(e[8]) for e in pulled.dict()),default=None)
    if kp<0:
        assert all(int(e[8])>=-kp for e in pulled.dict()),('WHOLE_RELATIVE_REGULARITY_FAILED',degree,rdegree,kp,min_kappa)
        regular=P({tuple(int(n)+(kp if j==8 else 0) for j,n in enumerate(e)):c for e,c in pulled.dict().items()})
    else:regular=kap**kp*pulled
    ans+=delta**dp*regular
    item={'ordinary_degree':int(degree),'R_degree':int(rdegree),'kappa_prefactor':int(kp),
          'delta_prefactor':int(dp),'whole_pullback_terms':before,
          'whole_pullback_min_kappa':min_kappa,'regular_terms':len(regular.dict())}
    records.append(item);layer_polys.append({'record':item,'whole_pullback':pulled,'regular':regular})
    print('layer',item,'seconds',time.monotonic()-st,flush=True)
    # Avoid reserializing the growing rational-function payload per layer.
print('COMPLETE G terms',len(ans.dict()),'seconds',time.monotonic()-st,flush=True)
monomial_decomposition=verify_and_export(P,ans,out)
# Every incidence coefficient denominator is checked at the protected
# cyclic ratio point without specializing either material parameter.
K0=FractionField(PolynomialRing(QQ,names=('kk','rr')));kk,rr=K0.gens()
sp=K.ring().hom([kk,rr,1,-1],K0);bad=[];denoms=set()
for c in ans.coefficients():
    den=sp(c.denominator())
    assert den!=0,'coefficient denominator vanishes at cyclic ratio section'
    for factor,power in K0(den).numerator().factor():
        factor=K0.ring()(factor)
        if factor.is_constant():continue
        monic=factor/factor.leading_coefficient();denoms.add(str(monic))
        assert any(monic==z for z in [K0.ring()(kk),K0.ring()(rr),K0.ring()(1+kk),K0.ring()(1+rr)])
# Compare the COMPLETE polynomial with the independently produced cyclic G.
cy=load(str(p/'independent_subcritical_single_RR_full_complete_E_chart_H1.sobj'))
CP=cy['ring'];CK=CP.base_ring();ck,cr=CK.gens();spc=K.ring().hom([ck,cr,1,-1],CK)
cy_image=CP({e:spc(c.numerator())/spc(c.denominator()) for e,c in ans.dict().items()})
target=next(v['source'] for v in cy['items'] if v['gauge']==2 and v['Q_index']==6)
assert cy_image==target,'complete cyclic polynomial mismatch'
report={'status':'PASS_COMPLETE_RELATIVE_REGULARITY','scope':__doc__,
        'parameters':['k','r','eta','zeta'],'source_hashes':{n:hashlib.sha256((p/n).read_bytes()).hexdigest() for n in [gn,en]},
        'gauge':gi,'Q_index':qi,'layers':records,'terms':len(ans.dict()),
        'all_negative_kappa_layers_exactly_divisible_before_any_quotient':True,
        'all_coefficient_denominators_protected_at_cyclic_section':True,
        'cyclic_denominator_material_factors':sorted(denoms),
        'complete_cyclic_polynomial_matches_independent_source':True,
        'complete_G_membership_in_t_kappa':monomial_decomposition,
        'seconds':time.monotonic()-st}
save({'ring':P,'source':ans,'layers':layer_polys,'report':report,'scope':__doc__},str(out/'generic_G_relative_regular.sobj'),compress=False)
(out/'generic_G_relative_regular_sparse.json').write_text(json.dumps({
    'coefficient_field':'QQ(k,r,eta,zeta)',
    'variables':list(P.variable_names()),
    'definition':'kappa^13*delta^6*G(p21,Q6) after the exact incidence substitution',
    'terms':[{'exponents':list(map(int,e)),'coefficient':str(c)} for e,c in ans.dict().items()]},indent=2)+'\n')
(out/'generic_G_relative_regular.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
