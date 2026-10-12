"""Independently reconstruct both cleared graph polynomials by six Taylor coefficients."""
from pathlib import Path
import time,json
from sage.all import load,save
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
d=load(str(out/'pencil_fast_numerators.sobj'))
old=load(str(base/'round9/rank1_ratio_characteristic_linear_pair.sobj'))
Flat=d['flat_ring'];CR=d['base_ring'];K=CR.fraction_field();a,q,T,Hb=map(K,CR.gens());r=q-1
phi=old['ring'].base_ring().hom([a,a*q/(1+a*r),q*T/r,Hb/r],K)
D=d['determinant'];Nx,Ny=d['cramer_numerators'];u,v,s=Flat.gens()[:3]
checks={}
def check(key,val):checks[key]=bool(val);assert checks[key],key
def coefficient(source,i,j,den):
    terms={}
    for mon,co in source.dict().items():
        if (mon[2],mon[4])!=(i,j):continue
        cc=phi(co);quot,rem=den.quo_rem(cc.denominator());assert not rem
        coeff=cc.numerator()*quot
        for bm,bc in coeff.dict().items():
            key=(mon[0],mon[1],mon[3])+tuple(bm)
            terms[key]=terms.get(key,0)+bc
    return Flat(terms)
for index,key in enumerate(('remaining_mixed_source','quartic_source_Q0')):
    source=old[key];record=d['records'][index]
    check(key+'_quadratic_in_xy',all(mon[2]+mon[4]<=2 for mon in source.dict()))
    den=d['source_denominators'][index]
    cc={ij:coefficient(source,*ij,den) for ij in ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))}
    rebuilt=D*D*cc[0,0]+D*(Nx*cc[1,0]+Ny*cc[0,1])+Nx*Nx*cc[2,0]+Nx*Ny*cc[1,1]+Ny*Ny*cc[0,2]
    check(key+'_exact_cleared_source_replay',rebuilt==d['rows'][index]*s**record['removed_s_power'])
    print(key,'PASS',round(time.monotonic()-start,2),flush=True)
report={'status':'PASS_ALL_PARAMETER_SOURCE_REPLAY','checks':checks,'elapsed_seconds':time.monotonic()-start,
        'scope':'Exact identity of both saved graph numerators with original H and Q0 sources. This supplies no coprimality, parameter-exception exclusion, or dimension conclusion.'}
(out/'pencil_graph_independent_reader.json').write_text(json.dumps(report,indent=2)+'\n')
save({'checks':checks,'source_file':d['source_file']},str(out/'pencil_graph_independent_reader.sobj'))
