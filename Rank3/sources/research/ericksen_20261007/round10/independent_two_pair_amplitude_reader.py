"""Independent twelve-axis classification, rebuilt from normalized sources."""
from pathlib import Path
from itertools import product
import json,time
from sage.all import load,save,matrix,prod
out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
norm=load(str(out/'general_curve_homogeneous_normalized.sobj'))
lead=load(str(out/'general_curve_dimensionless_leading.sobj'))
d=load(str(out/'general_curve_two_pair_offdiagonal.sobj'))
cert=load(str(out/'general_curve_two_pair_amplitude_certificate.sobj'))
R=d['ring'];K=R.base_ring();C=K.ring();H,Z,S,N=K.gens();x=R.gens()
guards=[C(H),C(Z),C(H-Z),C(Z-1),C(H+1),C(H-Z+1),C(S-1),C(H*S-Z+1),C((H+1)*S-Z)]
def isunit(c):
    if not c:return False
    for pol in (c.numerator(),c.denominator()):
        for fac,e in C(pol).factor():
            fac=C(fac)
            if fac.is_constant():continue
            if not any(all((K(fac)/g).derivative(v)==0 for v in C.gens()) for g in guards):return False
    return True
im=[x[0],-2*(H*x[1]+Z*x[0])/(1-S),x[2],-2*H*(x[3]+Z*x[2])/(1+H*S-Z),x[4],2*Z*(x[5]-H*x[4])/((1+H)*S-Z)]
checks['same_protected_coordinate_map']=im==d['images'] and all(isunit(c.denominator()) for p in im for c in p.coefficients())
ph=norm['ring'].hom(im,R)
checks['six_cubic_sources']=[ph(f) for f in norm['cubics']]==d['cubics']
checks['three_quartic_sources']=[ph(f) for f in norm['quartics']]==d['quartics']
combs=[(-1,1,0),(-1,1,0),(1,0,1)]
checks['pure_combinations']=combs==cert['pure_pair_combinations']
for i in range(3):
    phi=R.hom([x[j] if j//2==i else R.zero() for j in range(6)],R)
    checks['pure_pair_square_'+str(i)]=sum(k*phi(f) for k,f in zip(combs[i],d['quartics']))==(-4 if i in (0,2) else 4)*x[2*i]**2*x[2*i+1]**2
seen=[]
for layer in d['layers']:
    missing=layer['missing'];active=[i for i in range(3) if i!=missing]
    A=layer['ring'];xx=A.gens();ims=[A.zero()]*6
    for j,t in enumerate(xx):ims[2*active[j//2]+j%2]=t
    phi=R.hom(ims,A)
    cubics=[phi(f) for f in d['cubics'] if phi(f)]
    quartics=[phi(f) for f in d['quartics']]
    BB=matrix(K,2,2,[K(co(0,0,H,Z,S,N)) for co in lead['blocks'][tuple(active)].list()])
    mix=phi(sum(im[2*active[0]+i]*BB[i,j]*im[2*active[1]+j] for i in range(2) for j in range(2)))
    rows=[mix]+cubics+quartics
    checks['all_layer_sources_'+str(missing)]=rows==layer['polynomials']
    for i,j in product((0,1),repeat=2):
        key=str(missing)+'_'+str(i)+str(j);seen.append((missing,i,j))
        axis=A.hom([xx[k] if k in (i,2+j) else A.zero() for k in range(4)],A)
        polynomials=[axis(f) for f in rows]
        item=next(it for it in cert['certificates'] if it['missing']==missing and tuple(it['axes'])==(i,j))
        if 'source_position' in item:
            cc=polynomials[item['source_position']].monomial_coefficient(A.monomial(*tuple(item['monomial'])))
            checks['original_axis_coefficient_'+key]=cc==item['coefficient']
            checks['coefficient_unit_'+key]=isunit(cc)
        else:
            checks['only_two_exceptions_'+key]=(missing,i,j) in ((0,1,1),(1,1,0))
            wall=H*S-Z if missing==0 else S
            qco=polynomials[item['quartic_source_position']].monomial_coefficient(xx[i]**2*xx[2+j]**2)
            checks['quartic_wall_source_'+key]=qco==item['quartic_coefficient'] and wall==item['wall'] and isunit(qco/wall**2)
            sv=Z/H if missing==0 else K.zero()
            def spec(f,nn=N):return A({m:K(cc.numerator()(H,Z,sv,nn))/K(cc.denominator()(H,Z,sv,nn)) for m,cc in f.dict().items()})
            mixed=spec(polynomials[0]);co=mixed.monomial_coefficient(xx[i]*xx[2+j])
            checks['mixed_curvature_source_'+key]=mixed==item['mixed_after_wall'] and co/N==item['mixed_N_unit'] and isunit(co/N)
            checks['sufficiency_all_selected_sources_'+key]=all(not spec(f,K.zero()) for f in polynomials)
checks['complete_twelve_axes']=len(seen)==12 and len(cert['certificates'])==12
assert all(checks.values()),[k for k,v in checks.items() if not v]
save({'checks':checks,'scope':'Complete fixed-direction two-pair amplitude classification on all nine guards, then combined with the separately audited 145-check entire-family boundary cuts. Does not control varying directions.'},str(out/'independent_two_pair_amplitude_reader.sobj'))
report={'status':'EXACT_ALL_TWELVE_AMPLITUDE_AXIS_CASES_PASS','checks':checks,'check_count':len(checks),'elapsed_seconds':time.monotonic()-start}
(out/'independent_two_pair_amplitude_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='checks'}),flush=True)
