"""Independent original-minor and ordered-material audit of two hyperbolas."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent; started=time.monotonic(); checks=[]
def ck(name,value):
    if not bool(value): raise AssertionError(name)
    checks.append(name)
charts=load(str(p/'independent_ratio_canonical_partial_zero_hyperbolas.sobj'))
witness=load(str(p/'independent_ratio_partial_zero_unit_witness.sobj'))
res=load(str(p/'independent_ratio_partial_zero_residuals.sobj'))
old=load(str(p/'independent_ratio_canonical_partial_zero_minor_span.sobj'))
signs=load(str(p/'independent_ratio_partial_zero_factor_signs.sobj'))
K=charts['ring'].base_ring(); a,r=charts['parameters']; q=1+r; m=1+a*r
P=K.ring(); U=PolynomialRing(K,names=('audit_x',)); x=U.gen()
V=PolynomialRing(K,names=('audit_z',)); z=V.gen()
def cx(f): return U({(int(e[0]),):K(c) for e,c in f.dict().items()})
def cz(f): return V({(int(e[0]),):K(c) for e,c in f.dict().items()})
def determinant_dp(M):
    # Independently expand the exterior product, without division or pivots.
    n=M.nrows(); d={0:U.one()}
    for i in range(n):
        nd={}
        for mask,c in d.items():
            for j in range(n):
                if (mask>>j)&1 or not M[i,j]: continue
                sg=(-1)**sum((mask>>k)&1 for k in range(j+1,n))
                key=mask|(1<<j); nd[key]=nd.get(key,U.zero())+sg*c*M[i,j]
        d=nd
    return d.get((1<<n)-1,U.zero())
def strip(f):
    if not f:return V.zero(),0
    e0=min(int(e) for e in f.dict())
    ck('strip_even_'+str(len(checks)),all((int(e)-e0)%2==0 for e in f.dict()))
    return V({(int(e)-e0)//2:c for e,c in f.dict().items()}),e0
T=PolynomialRing(QQ,names=('gap_u','gap_v')); u,v=T.gens(); TF=T.fraction_field()
ap=TF((u+v)/u); rp=TF(-u/(1+u+v))
def positive_pull(f):
    return sum(TF(c)*ap**int(e[0])*rp**int(e[1]) for e,c in P(f).dict().items())
ck('ordered_inverse_gap_u',(-r/m)/(1+a*(-r/m))==-r)
ck('ordered_inverse_a',((-r/m)+(a-1)*(-r/m))/(-r/m)==a)
records=[]
for wi in witness['items']:
    case=int(wi['case']); ci=charts['items'][case]; ri=res['items'][case]; oi=old['items'][case]
    M=matrix(U,[[cx(f) for f in row] for row in ci['matrix'].rows()])
    om=int(wi['old_minor_omitted_row']); ids=[j for j in range(7) if j!=om]
    genuine=determinant_dp(M.matrix_from_rows(ids)); g,ep=strip(genuine)
    ck('original_old6minor_'+str(case),genuine==cx(oi['original_minors'][om]))
    ck('old6minor_stripping_'+str(case),ep==int(wi['old_minor_divided_x_power']))
    ck('old6minor_polynomial_'+str(case),g==cz(wi['polynomials'][0]))
    seed=list(map(int,ri['seed'])); B=M.matrix_from_rows(seed)
    # The alternating cofactor vector is independently formed from all six
    # maximal minors of B. Row expansion gives det([B;row])=row.wedge.
    wedge=vector(U,[(-1)**(5+j)*determinant_dp(B.matrix_from_columns([k for k in range(6) if k!=j])) for j in range(6)])
    for j in range(6):ck('cramer_wedge_%d_%d'%(case,j),wedge[j]==cx(ri['cramer_numerators'][j]))
    ck('wedge_annihilates_seed_'+str(case),all(not f for f in B*wedge))
    polys=[g]; powers=[ep]
    for k,rn in enumerate(map(int,wi['cramer_residual_source_indices'])):
        raw=determinant_dp(matrix(U,list(B.rows())+[M.row(rn)])); pp,ee=strip(raw)
        ck('true_appended6minor_%d_%d'%(case,rn),raw==cx(ri['original_residuals'][rn]))
        ck('residual_cofactor_%d_%d'%(case,rn),raw==M.row(rn)*wedge)
        ck('residual_strip_power_%d_%d'%(case,rn),ee==int(ri['stripped_x_powers'][rn]))
        ck('residual_polynomial_%d_%d'%(case,rn),pp==cz(wi['polynomials'][k+1]))
        polys.append(pp);powers.append(ee)
    cofs=[cz(f) for f in wi['coefficient_polynomials']]
    ck('original_three_source_identity_'+str(case),sum(c*f for c,f in zip(cofs,polys))==1)
    den=P(wi['common_denominator']); actual=lcm([K(c).denominator() for f in cofs for c in f.coefficients()])
    ck('common_denominator_'+str(case),den==actual)
    integer_cofs=[]
    for j,f in enumerate(cofs):
        for c in f.coefficients():ck('cofactor_clear_%d_%d_%d'%(case,j,len(checks)),K(den*c).denominator()==1)
        integer_cofs.append(den*f)
    ck('cleared_original_identity_'+str(case),sum(c*f for c,f in zip(integer_cofs,polys))==den)
    factors=[]
    for fac,power in den.factor():
        fac=P(fac); pulled=positive_pull(fac); ns={sign(c) for c in pulled.numerator().coefficients()}; ds={sign(c) for c in pulled.denominator().coefficients()}
        same=(len(ns)==1 and len(ds)==1)
        # Ordered gap variables cover a>1, -1/a<r<0. Every denominator
        # factor has a numerator and denominator of constant coefficient sign.
        ck('ordered_factor_nonzero_%d_%d'%(case,len(factors)),same)
        if fac.total_degree()>4:
            linked=[]
            for rec in signs['items'][case]['tests']:
                rat=K(fac)/K(rec['factor'])
                if rat.numerator().total_degree()==0 and rat.denominator().total_degree()==0:
                    linked.append(rec)
                    val=positive_pull(rec['factor'])
                    ck('saved_gap_numerator_%d_%d'%(case,len(linked)),val.numerator()==T(rec['transformed']['numerator']))
                    ck('saved_gap_denominator_%d_%d'%(case,len(linked)),val.denominator()==T(rec['transformed']['denominator']))
            ck('largefactor_link_'+str(case),bool(linked))
        factors.append({'factor':str(fac),'power':int(power),'degree':int(fac.total_degree()),'gap_numerator_terms':len(pulled.numerator().dict()),'numerator_sign':int(next(iter(ns))),'denominator_sign':int(next(iter(ds)))})
    records.append({'case':case,'old_omitted':om,'seed':seed,'residual_rows':list(wi['cramer_residual_source_indices']),'stripped_x_powers':powers,'denominator':str(den.factor()),'factors':factors,'cofactor_degrees':[int(f.degree()) for f in cofs]})
    print('CASE',case,'ORIGINAL SIX-MINORS / UNIT / ALL ORDERED FACTORS PASS',len(checks),round(time.monotonic()-started,2),flush=True)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-started,'records':records,'scope':'Cases0/1 full normal rank6 at every complex point of their off-axis hyperbolas for every ordered positive simple spectrum. Exact combinations use only actual six-minors and x units, not a non-tau pivot or the normalized gcd old_quad. All denominator factors have independently proved strict sign via positive ordered gap variables. Case2 not included.'}
save({'checks':checks,'records':records},str(p/'independent_partial_zero_unit_reader.sobj'))
(p/'independent_partial_zero_unit_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(report,flush=True)
