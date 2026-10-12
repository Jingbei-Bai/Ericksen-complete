"""Whole-parameter Cramer substitution using polynomial arithmetic over QQ.

The common p23 factor in det, Nx and Ny is removed before substitution.
No repeated arithmetic in a rational function coefficient field is used
inside the graph substitution. No determinant norm is discarded.
"""
from pathlib import Path
import time,json,os
from sage.all import QQ,PolynomialRing,load,save
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
d=load(str(base/'round9/rank1_ratio_characteristic_linear_pair.sobj'))
oldR=d['ring'];oldK=oldR.base_ring()
CR=PolynomialRing(QQ,names=('a','q','T','Hb'));K=CR.fraction_field();a,q,T,Hb=map(K,CR.gens())
r=q-1;b=a*q/(1+a*r);c=q*T/r;h=Hb/r
phi=oldK.hom([a,b,c,h],K)
R=PolynomialRing(K,names=('u','v','s'));Flat=PolynomialRing(QQ,names=('u','v','s','a','q','T','Hb'))
def transform(p):return {tuple(m):phi(co) for m,co in p.dict().items()}
raw=[transform(d['determinant'])]+[transform(p) for p in d['cramer_numerators']]
basepolys=[]
for row in raw:
    assert all(m[2]==m[4]==0 and m[3]>=1 for m in row)
    basepolys.append({(m[0],m[1],m[3]-1):co for m,co in row.items()})
assert d['determinant_s_power']==1
def lcm_den(rows):
    result=CR.one()
    for row in rows:
        for co in row.values():result=result.lcm(co.denominator());result=CR(result/result.lc())
    return result
def clear_coeff(co,den):
    quo,rem=den.quo_rem(co.denominator());assert not rem
    return co.numerator()*quo
def flatten(row,den):
    ans={}
    for m,co in row.items():
        for pe,cc in clear_coeff(co,den).dict().items():ans[tuple(m)+tuple(pe)]=cc
    return Flat(ans)
GD=lcm_den(basepolys);D,Nx,Ny=[flatten(row,GD) for row in basepolys]
print('prepared',round(time.monotonic()-start,2),'flat_terms',[p.number_of_terms() for p in (D,Nx,Ny)],flush=True)
up,vp,sp=Flat.gens()[:3]
def base_flat(p):return Flat({(0,0,0)+tuple(m):co for m,co in p.dict().items()})
sources=[('H',d['remaining_mixed_source']),('Q0',d['quartic_source_Q0'])]
rows=[];records=[];denominators=[]
for name,source in sources:
    row=transform(source);GS=lcm_den([row]);power=max(m[2]+m[4] for m in row)
    powersD=[D**j for j in range(power+1)];powersX=[Nx**j for j in range(power+1)];powersY=[Ny**j for j in range(power+1)]
    ans=Flat.zero()
    for m,co in row.items():
        mon=up**m[0]*vp**m[1]*sp**m[3]
        ans+=base_flat(clear_coeff(co,GS))*mon*powersX[m[2]]*powersY[m[4]]*powersD[power-m[2]-m[4]]
    s_power=min(m[2] for m in ans.dict()) if ans else 0
    if s_power:ans=Flat({tuple(e-(s_power if j==2 else 0) for j,e in enumerate(m)):co for m,co in ans.dict().items()})
    parity=sorted(set(tuple(e%2 for e in mon[:3]) for mon in ans.dict()))
    rows.append(ans);denominators.append(GS)
    rec={'source':name,'det_clear_power':int(power),'removed_s_power':int(s_power),
         'flat_terms':ans.number_of_terms(),'joint_degree':max(sum(m[:3]) for m in ans.dict()),
         'joint_support':len(set(tuple(m[:3]) for m in ans.dict())),'parities':parity}
    records.append(rec);print(json.dumps(rec),'elapsed',round(time.monotonic()-start,2),flush=True)
    payload={'flat_ring':Flat,'base_ring':CR,'fiber_ring':R,'parameters':(a,q,T,Hb),
       'parameter_map':{'b':b,'c':c,'h':h},'determinant':D,'cramer_numerators':[Nx,Ny],
       'common_graph_denominator':GD,'source_denominators':denominators,'rows':rows,'records':records,
       'source_file':'round9/rank1_ratio_characteristic_linear_pair.sobj',
       'identity':'row_k * s^removed = GS_k * D^power * source_k(x=Nx/D,y=Ny/D), where D,Nx,Ny share the same cleared coefficient denominator; graph variables are u,v,s.',
       'scope':'Exact graph numerators on original p23*det!=0 chart. No gcd, coprimality, dimension, or spectrum-exception claim.'}
    tmp=out/'pencil_fast_numerators_writing.sobj';save(payload,str(tmp));os.replace(tmp,out/'pencil_fast_numerators.sobj')
    (out/'pencil_fast_numerators.json').write_text(json.dumps({'status':'PARTIAL' if len(rows)<2 else 'COMPLETE_EXACT_NUMERATORS','records':records,'elapsed_seconds':time.monotonic()-start,'scope':payload['scope']},indent=2)+'\n')
print('DONE',round(time.monotonic()-start,2),flush=True)
