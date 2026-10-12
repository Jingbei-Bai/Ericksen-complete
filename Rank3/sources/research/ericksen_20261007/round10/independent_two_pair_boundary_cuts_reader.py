"""Independent full-source replay of both uniform two-pair amplitude cuts."""
from pathlib import Path
import time,json
from sage.all import load,save,QQ,PolynomialRing,matrix,vector,prod

out=Path(__file__).resolve().parent;start=time.monotonic();checks={};cert=[]
full=load(str(out/'general_curve_dimensionless_full.sobj'))
extra=load(str(out/'general_curve_extra_quartic_cubics.sobj'))
quart=load(str(out/'general_curve_remaining_quartics.sobj'))
original=full['rows']+extra['rows']+quart['rows']

def units(poly,guards):
    C=guards[0].parent();poly=C(poly);K=C.fraction_field()
    if not poly:return False
    for f,e in poly.factor():
        f=C(f)
        if f.is_constant():continue
        if not any(all((K(f)/g).derivative(v)==0 for v in C.gens()) for g in guards):return False
    return True

for proportional in [False,True]:
    stem='general_curve_proportional_Z' if proportional else 'general_curve_constant_Z'
    d=load(str(out/(stem+'_full_sources.sobj')))
    kr=load(str(out/(stem+'_quadratic_relations.sobj')))
    cut=load(str(out/(stem+'_boundary_cut.sobj')))
    R=d['ring'];K=R.base_ring();C=K.ring();a,r,H,c=K.gens();xs=R.gens()
    Z=c*H if proportional else c;S=c if proportional else K.zero();args=(a,r,H,Z,S,0)
    guards=[C(a),C(r),C(1+r),C(1+a*r),C(H),C(Z),C(Z-H),C(1-Z),C(1+H),C(1+H-Z),C(1-S),C(1+H*S-Z),C((1+H)*S-Z)]
    images=[xs[0],-2*(H*xs[1]+Z*xs[0])/(1-S),xs[2],-2*H*(xs[3]+Z*xs[2])/(1+H*S-Z),xs[4],2*Z*(xs[5]-H*xs[4])/((1+H)*S-Z)]
    label='proportional' if proportional else 'constant'
    checks[label+'_images']=images==d['images'] and len(d['rows'])==16
    # Verify the image directly inverts the original opposite-pair graph.
    w1,m1,w2,m2,w3,m3=images
    checks[label+'_original_opposite_pairs']=[w1,-((1-S)*m1/2+Z*w1)/H,w2,-(1+H*S-Z)*m2/(2*H)-Z*w2,w3,((1+H)*S-Z)*m3/(2*Z)+H*w3]==list(xs)
    coefficients={};monomials={}
    def coeff(cc):
        if cc not in coefficients:coefficients[cc]=K(cc.numerator()(*args))/K(cc.denominator()(*args))
        return coefficients[cc]
    def pull(p):
        ans=R.zero()
        for mon,cc in p.dict().items():
            mon=tuple(mon)
            if mon not in monomials:monomials[mon]=prod(images[i]**mon[i] for i in range(6))
            ans+=coeff(cc)*monomials[mon]
        return ans
    for i,p in enumerate(original):
        checks[label+'_original_full_source_'+str(i)]=pull(p)==d['rows'][i]
        con,den=d['scalings'][i]
        checks[label+'_primitive_replay_'+str(i)]=d['primitive_rows'][i]*K(con)/den==d['rows'][i]
        checks[label+'_source_denominators_'+str(i)]=all(units(c.denominator(),guards) for c in d['rows'][i].coefficients())
    print(label,'16 full sources PASS',round(time.monotonic()-start,2),flush=True)
    polys=[]
    for pos,pair in enumerate(((0,1),(0,2),(1,2))):
        other=3-sum(pair)
        polys.extend(d['rows'][pos]*xs[j] for j in [2*other,2*other+1])
    polys.extend(d['rows'][3:7])
    monomials3=[xs[i]*xs[j+2]*xs[k+4] for i in range(2) for j in range(2) for k in range(2)]
    mat=matrix(K,[[f.monomial_coefficient(mon) for f in polys] for mon in monomials3])
    checks[label+'_full_ten_sources']=polys==kr['source_polys']
    checks[label+'_symbol_matrix']=mat==kr['tensor_matrix']
    checks[label+'_generic_kernel_dimension']=mat.rank()==6 and kr['kernel'].nrows()==4
    multipliers=[]
    for i,(v,den) in enumerate(zip(kr['kernel'].rows(),kr['kernel_denominators'])):
        cs=[K(den)*vv for vv in v];f=sum((cc*p for cc,p in zip(cs,polys)),R.zero())
        checks[label+'_cleared_multipliers_'+str(i)]=all(cc.denominator().is_constant() for cc in cs)
        checks[label+'_symbol_identity_'+str(i)]=not mat*vector(cs)
        checks[label+'_full_combination_'+str(i)]=f==kr['quadratics'][i] and f.total_degree()<=2
        multipliers.append(cs)
    f=kr['quadratics'][0];checks[label+'_cut_source']=f==cut['full_quadratic']
    kept=(3,5) if proportional else (1,4)
    restriction=R({mon:co for mon,co in f.dict().items() if sum(mon)==2 and all(mon[j]==0 for j in range(6) if j not in kept)})
    co=f.monomial_coefficient(xs[kept[0]]**2);ratio=(1+r)*c**2 if proportional else (1+r)*c**2/r
    checks[label+'_line_restriction']=restriction==cut['bad_line_restriction']==co*(xs[kept[0]]**2+ratio*xs[kept[1]]**2)
    unit=4*c*(1+r)*(c-1)**2*H**2/r if proportional else c*(1+r)*(c-1)**2/2
    pp=C(co/unit)
    checks[label+'_coefficient_identity']=unit==cut['protected_unit'] and co==cut['factor_coefficient'] and pp==cut['parameter_polynomial']
    checks[label+'_unit_protected']=units(unit.numerator(),guards) and units(unit.denominator(),guards)
    checks[label+'_independent_spectra']=not pp.derivative(C.gen(0)) and not pp.derivative(C.gen(1))
    if proportional:
        checks[label+'_every_curve_parameter']=pp.subs({C.gen(2):0})==32*(1-C.gen(3))
        checks[label+'_cminusone_survives']=bool(pp.subs({C.gen(3):-1})) and all(not cc.denominator().subs({C.gen(3):-1})==0 for cc in multipliers[0])
        # Original source denominators remain protected and nonzero as rational
        # functions of H at c=-1; the cleared multipliers have no c+1 inverse.
        checks[label+'_full_identity_at_cminusone']=sum((cc.subs({C.gen(3):-1})*p.map_coefficients(lambda x:x.subs({C.gen(3):-1})) for cc,p in zip(multipliers[0],polys)),R.zero())==f.map_coefficients(lambda x:x.subs({C.gen(3):-1}))
    else:
        U=PolynomialRing(QQ,'z');zv=U.gen();cs=[U.zero()]*(pp.degree(C.gen(2))+1)
        for mon,cc in pp.dict().items():
            assert not mon[0] and not mon[1]
            cs[mon[2]]+=cc*zv**mon[3]
        checks[label+'_H_coefficients']=cs==cut['H_coefficients']
        gg=sum(mm*cc for mm,cc in zip(cut['Bezout_multipliers'],cs))
        checks[label+'_Bezout_source_identity']=gg==cut['coefficient_gcd']==1
    cert.append({'family':label,'multipliers':multipliers,'symbol_matrix':mat,'first_quadratic':f,'restriction':restriction,'parameter_polynomial':pp})
    print(label,'full combinations and cut PASS',round(time.monotonic()-start,2),flush=True)
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'certificates':cert,'scope':'Uniform cuts of the two stated amplitude P1 boundaries by full necessary sources, for generic H at each permitted spectrum and each permitted fixed curve constant. Cleared identities include kernel pivot walls. Does not by itself classify every projective boundary.'},str(out/'independent_two_pair_boundary_cuts_reader.sobj'))
report={'status':'EXACT_TWO_ENTIRE_CURVE_FAMILY_BOUNDARY_CUTS_PASS','check_count':len(checks),'checks':checks,'elapsed_seconds':time.monotonic()-start}
(out/'independent_two_pair_boundary_cuts_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
