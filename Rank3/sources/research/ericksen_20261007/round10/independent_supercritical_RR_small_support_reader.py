"""Independent full-field reconstruction of the no-D supercritical certificate.

Unlike the creator, form the differential operator from joint_fields_direct
minus the full Q0 gauge; differentiate the six original full R sources.
Degree-support enumeration proves that only R12 can contribute to H13.
Read the stored signed DAG only for the already audited H32 module,
never to obtain the new H13 comparison.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;start=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
s=load(str(p/'independent_supercritical_RR_small_support_sources.sobj'))
g=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
fd=load(str(p/'independent_ratio_weighted_degreefall.sobj'))
raw=load(str(p/'joint_fields_direct.sobj'))
red=load(str(p/'independent_ratio_diagonal_Q0_reduction.sobj'))
mod=load(str(p/'independent_delta_H_global_L_module.sobj'))
bd=load(str(p/'independent_ratio_RR_boundary_leading.sobj'))
cl=load(str(p/'independent_ratio_RR_boundary_classification.sobj'))
ex=load(str(p/'independent_ratio_two_RR_mixed_exclusion.sobj'))
old=raw['ring'].base_ring();a0,b0,h0,z0=old.gens()
P=PolynomialRing(old,names=tuple('chk_p%d'%i for i in range(6))+tuple('chk_R%d'%i for i in range(3)))
pv=P.gens();pool=g['polynomial_pool'];wt=(1,)*6+(2,)*3
def cp(f):return P({tuple(map(int,e)):old(c) for e,c in f.dict().items() if c})
def pp(dd):return P({tuple(map(int,e)):old(c) for e,c in dd.items() if c})
def ref(j):return pp(pool[j])
def bi(e):return sum(e[:6]),sum(e[6:])
def layer(f,ij):return P({tuple(map(int,e)):c for e,c in f.dict().items() if bi(tuple(map(int,e)))==ij})
dcache={}
def derivcoeff(f,par):
    out={}
    for e,c in f.dict().items():
        key=(c,par)
        if key not in dcache:dcache[key]=c.derivative(par)
        if dcache[key]:out[tuple(map(int,e))]=dcache[key]
    return P(out)
def productface_all(pairs,ij):
    ans={}
    for scalar,f1,f2 in pairs:
        for ee,c in f1.dict().items():
            e=tuple(map(int,ee));be=bi(e)
            for ff,d in f2.dict().items():
                f=tuple(map(int,ff));bf=bi(f)
                if (be[0]+bf[0],be[1]+bf[1])!=ij:continue
                ne=tuple(e[j]+f[j] for j in range(9));ans[ne]=ans.get(ne,old.zero())+scalar*c*d
    return P({e:c for e,c in ans.items() if c})
R=[cp(it['full']) for it in fd['items']];Q0=cp(fd['Q0'])
F=[];base=[]
for i in range(3):
    ff=[cp(raw['fields'][i,j]) for j in range(9)];ff[6+i]-=old(red['items'][i]['kappa'])*Q0
    F.append(ff);base.append([cp(raw['base_fields'][i,j]) for j in range(2)])
    for j in range(9):ck('actual field full replay %d %d'%(i,j),ff[j]==ref(g['full_first_gauged_field_refs'][i][j]))
    for j in range(2):ck('actual base full replay %d %d'%(i,j),base[i][j]==ref(g['full_base_field_refs'][i][j]))
    ck('full field bounds '+str(i),all(all(sum(e)<= (3 if j==6+i else 2) and sum(w*n for w,n in zip(wt,e))<=wt[j]+2 for e in f.dict()) for j,f in enumerate(ff)))
    ck('base bounds '+str(i),all(all(sum(e)<=1 for e in f.dict()) for f in base[i]))
D=matrix(P,[[ref(j) for j in row] for row in g['D_refs']]);C=matrix(P,[[ref(j) for j in row] for row in g['C_refs']]);Lmat=matrix(P,[[ref(j) for j in row] for row in g['L_refs']]);ell=Lmat*vector(P,pv[6:])
ZZ=vector(P,[ref(j) for j in g['Z_source_refs']]);DZ=D*ZZ;DCe=D*C*ell
for j,f in enumerate(R):
    ck('full original R '+str(j),f==ref(g['full_source_refs'][j]))
    ck('R ordinary weighted bound '+str(j),all(sum(e)<=4 and sum(w*n for w,n in zip(wt,e))<=5 for e in f.dict()))
    ck('R12 complete source '+str(j),layer(f,(1,2))==pp(s['R12_original'][j]))
    ck('R31 global ell module '+str(j),layer(f,(3,1))==DCe[j])
    ck('Rtop factor identity '+str(j),ref(g['source_top4_refs'][j])==DZ[j])
H13=[]
for k,rec in enumerate(g['records']):
    i=int(rec['direction']);j=int(rec['source_position']);f=R[j]
    correct=[ref(n) for n in rec['correction_coefficients']]
    # Enumerate ALL full supports before selecting the small derivative.
    possible=True
    for ef in f.dict():
        ef=tuple(map(int,ef))
        for l in range(9):
            if not ef[l]:continue
            for ev in F[i][l].dict():
                ee=tuple(ef[n]+int(ev[n])-(1 if n==l else 0) for n in range(9))
                if bi(ee)==(1,3):possible=possible and bi(ef)==(1,2)
        for bf in base[i]:
            for ev in bf.dict():
                if bi(tuple(ef[n]+int(ev[n]) for n in range(9)))==(1,3):possible=possible and bi(ef)==(1,2)
    ck('full support only R12 contributes H13 '+str(k),possible)
    ck('full polynomial correction has no H13 '+str(k),all((sum(e[:6])+sum(fx[:6]),sum(e[6:])+sum(fx[6:]))!=(1,3) for l in range(6) for e in correct[l].dict() for fx in R[l].dict()))
    f12=layer(f,(1,2));df=[f12.derivative(v) for v in pv];db=[derivcoeff(f12,t) for t in (h0,z0)]
    pairs=[(old.one(),F[i][l],df[l]) for l in range(9)]
    pairs += [(old.one(),base[i][l],db[l]) for l in range(2)]
    hh=productface_all(pairs,(1,3));H13.append(hh)
    ck('full differential H13 '+str(k),hh==pp(s['H13_original'][k]))
    # Undo the universal binomial normalization coefficient by coefficient.
    mi=mod['items'][k];sc=mod['original_to_binomial_variable_scales'];co=[]
    for poly in mi['cofactors']:
        co.append(P({tuple(map(int,e)):old(c)/prod(sc[l]**int(e[l]) for l in range(9)) for e,c in poly.dict().items()}))
    # This full-DAG module was independently linked to the genuine source
    # by independent_intermediate_atlas_H_reader (932 checks). Replay its
    # original-coordinate polynomial identity here without redifferentiating
    # irrelevant coefficients of the complete R a second time.
    norm=mi['normalized_face'];sc=mod['original_to_binomial_variable_scales']
    h32=P({tuple(map(int,e)):old(c)/prod(sc[l]**int(e[l]) for l in range(9)) for e,c in norm.dict().items()})
    ck('H32 complete global ell identity '+str(k),h32==sum((co[l]*ell[l] for l in range(3)),P.zero()))
    ck('H32 cofactor bidegrees '+str(k),all(all(bi(tuple(map(int,e)))==(2,1) for e in c.dict()) for c in co))
    # Rtop=DZ provides the only possible ordinary-degree-six cancellation.
    top=ref(g['source_top4_refs'][j]);zsrc=ref(g['Z_source_refs'][int(g['sigma'][i])])
    sig=int(g['sigma'][i]);cc=vector(P,[sum((correct[l]*D[l,jj] for l in range(6)),P.zero()) for jj in range(3)])
    ck('all degree-six terms cancel by exact source factors '+str(k),cc[sig]==top.derivative(pv[6+i]) and all(cc[jj]==0 for jj in range(3) if jj!=sig))
    print('independent global H',k,'PASS',round(time.monotonic()-start,2),flush=True)
K=bd['coefficient_field'];a,r,H,Z=bd['parameters'];q=1+r;G=bd['physical_ring'];ph=old.ring().hom(bd['parameter_map'],K);cache={}
def cv(c):
    c=old(c)
    if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
    return cache[c]
def mapped(f):return G({tuple(map(int,e[6:]))+tuple(map(int,e[:6])):cv(c) for e,c in f.dict().items()})
for m,it in enumerate(s['items']):
    ch=it['chart'];T=it['H13_values'][0].parent();v,x,y=T.gens();ev=G.hom(list(ch['RR'])+list(ch['p']),T);sq=ch['v_squared']
    def rv(f):
        dd={}
        for ee,c in T(f).dict().items():
            e=tuple(map(int,ee));ne=(e[0]%2,e[1],e[2]);dd[ne]=dd.get(ne,K.zero())+c*sq**(e[0]//2)
        return T({e:c for e,c in dd.items() if c})
    actual=[rv(ev(mapped(f))) for f in H13]
    for j,val in enumerate(actual):ck('support H13 independent pull %d %d'%(m,j),val==it['H13_values'][j])
    for j,cert in enumerate(it['amplitude_certificates']):
        pos=list(g['source_order']).index(int(cert['source_name'][1:]));picked=[actual[i*6+pos] for i in cert['directions']]
        ck('all endpoint amplitude Bezout %d %d'%(m,j),sum((a*b for a,b in zip(picked,cert['bezout_coefficients'])),T.fraction_field().zero())==cert['amplitude'])
    # No p amplitude is inverted in the units of these identities.
    for j,cert in enumerate(it['amplitude_certificates']):
        for unit in cert['source_units']:
            for po in (unit.numerator(),unit.denominator()):
                ck('unit monomial and no amplitude '+str(len(checks)),len(po.dict())==1 and all(tuple(map(int,e[1:]))==(0,0) for e in po.dict()))
                for c in po.coefficients():
                    for pol in (c.numerator(),c.denominator()):
                        for fac,mul in pol.factor():ck('unit base protected '+str(len(checks)),ex['items'][m]['guard'].quo_rem(fac)[1]==0)
    a1,a2=ch['necessary_nonzero_factors'];ck('two RR chart pivot cover '+str(m),2*a1-a2==(-3*H,-3,3)[m])
# Complete single-RR kernels and all zero-p endpoints: extract R12 only.
T=PolynomialRing(K,names=('ax0','ax1','ax2','ax3'));x,y,z,w=T.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));GG=H*Z-H+Z
guard=K.ring()(a*r*q*(a-1)*(1+a*r)*H*Z*(H-Z)*(H+1)*(Z-1)*(1+H-Z)*GG)
axes=[((1,0,0),(0,0,x,y,z,w),{10:z*(H+2),11:w*(H+2),13:x*(Z-2),19:y*(Z-2)}),((0,1,0),(x,y,0,0,z,w),{10:z*(2*H+1),11:w*(2*H+1),15:x*(Z-2*H),18:y*(Z-2*H)}),((0,0,1),(x,y,z,w,0,0),{13:z*(2*Z-1),15:x*(2*Z-H),18:y*(2*Z-H),19:w*(2*Z-1)})]
for m,(rr,ps,expected) in enumerate(axes):
    ev=G.hom(list(rr)+list(ps),T);Lglobal=bd['linear_forms'];ck('single RR ell zero '+str(m),all(ev(f)==0 for f in Lglobal))
    J=matrix(T,[[ev(f.derivative(vv)) for vv in G.gens()[3:]] for f in Lglobal]);B=matrix(K,[[T(f).monomial_coefficient(vv) for vv in T.gens()] for f in ps])
    ck('complete single RR kernel '+str(m),J.rank()==2 and B.rank()==4)
    for j,idx in enumerate(g['source_order']):
        val=ev(mapped(layer(R[j],(1,2))))
        if idx not in expected:ck('zero single R12 %d %d'%(m,idx),val==0);continue
        term=T(expected[idx]);var=next(vv for vv in T.gens() if term.degree(vv)==1);unit=K(val.monomial_coefficient(var))/K(term.monomial_coefficient(var))
        ck('single source complete factor %d %d'%(m,idx),unit!=0 and val==unit*term)
        for pol in (unit.numerator(),unit.denominator()):
            for fac,mul in pol.factor():ck('single source guard '+str(len(checks)),guard.quo_rem(fac)[1]==0)
walls=(H+2)*(Z-2)*(2*H+1)*(Z-2*H)*(2*Z-1)*(2*Z-H);sub=K.ring().hom([a,r,1,-1],K)
ck('one base open all spectra',sub(K.ring()(walls))==729 and sub(K.ring()(GG))==-3)
report={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-start,'scope':'All full actual fields and sources, R31/H32 global ell identities, direct H13 differential reconstruction, all54 two-RR evaluations and amplitude identities, every single-RR kernel and R12 factor. No D-rank assumption or p denominator. Supercritical no-lift follows from finite-e boundedness and valuation, not an assumed invariant support.'}
save({'checks':checks,'report':report},str(p/'independent_supercritical_RR_small_support_reader.sobj'))
(p/'independent_supercritical_RR_small_support_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
