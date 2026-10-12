"""Independent direct Taylor extraction and Leibniz determinant replay."""
from pathlib import Path
from itertools import permutations
from math import comb
import json,time
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_collision_resultant_series.sobj'));original=load(str(out/'pencil_vzero_sone_cubics.sobj'))
S=d['ring'];K=S.fraction_field();q,T=map(K,S.gens());r=q-1;a=d['a_on_wall'];root=d['root']
checks={}
def check(key,val):checks[key]=bool(val);assert checks[key],key
def param_map(co):return K(co(a,q,T,0))
cluster=[]
for index,(poly,height,den) in enumerate(zip(original['cubics'],(20,23),d['local_clearing_denominators'])):
    vv=[]
    for ell in reversed(range(4)):
        vals=[]
        for n in range(4):
            val=K.zero()
            for j in range(ell,4):
                hb_power=height-j-2+ell-n
                cc=original['base_ring']({tuple(0 if k==3 else ee for k,ee in enumerate(mon)):num
                                         for mon,num in poly[j].dict().items() if mon[3]==hb_power})
                val+=comb(j,ell)*param_map(cc)*root**(j-ell)
            cleared=val*K(den);check(f'cleared_{index}_{ell}_{n}',cleared.denominator()==1)
            vals.append(S(cleared))
        vv.append(vals)
    check(f'direct_Taylor_source_{index}',vv==d['cluster_coefficients'][index]);cluster.append(vv)
guard=S(d['parameters'][0])*S(d['parameters'][1])*(S(d['parameters'][0])-1)
for i,den in enumerate(d['local_clearing_denominators']):
    check(f'clearing_denominator_{i}_only_q_T_r',all(not guard.quo_rem(f)[1] for f,e in den.factor()))
print('direct cubic to cluster replay PASS',round(time.monotonic()-start,2),flush=True)
zero=[S.zero()]*4;matrix=[]
for vec in cluster:
    for shift in range(3):matrix.append([zero]*shift+vec+[zero]*(2-shift))
def mul(x,y):return [sum((x[j]*y[k-j] for j in range(k+1)),S.zero()) for k in range(4)]
coeff=[S.zero()]*4;count=0
for perm in permutations(range(6)):
    if any(not any(matrix[i][perm[i]]) for i in range(6)):continue
    value=[S.one(),S.zero(),S.zero(),S.zero()]
    for i in range(6):value=mul(value,matrix[i][perm[i]])
    sign=(-1)**sum(perm[i]>perm[j] for i in range(6) for j in range(i+1,6))
    for k in range(4):coeff[k]+=sign*value[k]
    count+=1
check('full_Leibniz_vs_saved_three_coefficients',coeff==d['coefficients'])
local=load(str(out/'pencil_collision_independent_reader.sobj'))
den_product=K(d['local_clearing_denominators'][0])**3*K(d['local_clearing_denominators'][1])**3
check('t4_global_local_prefactor',K(coeff[1])/den_product==local['global_t4'])
source=load(str(out/'pencil_collision_coefficient_ideal.sobj'))
qq,tt=S.gens();g=qq*tt*(tt-1)*(tt-2)*(qq-1)*((qq-1)*tt-4)
for i,co in enumerate(coeff[1:]):
    factor=co.factor();keep=S.one();removed=S.one()
    for ff,power in factor:
        if not g.quo_rem(ff)[1]:removed*=ff**power
        else:keep*=ff**power
    check(f'coefficient_{i}_strip_protected_only',keep==source['source_polynomials'][i] and keep*removed*factor.unit()==co)
nextdata=load(str(out/'pencil_collision_next_remainder.sobj'))
nf=[p for p,e in nextdata['next_obstruction'].numerator().factor() if p.degree(qq)>0];assert len(nf)==1
# The independent remainder obstruction is only compared ON P=0.
# It need not literally equal the full coefficient, which also contains
# terms proportional to the preceding resultant coefficient.
literal_remainder_equality=(source['source_polynomials'][1]/source['source_polynomials'][1].lc()==nf[0]/nf[0].lc())
report={'status':'PASS_INDEPENDENT_FULL_COLLISION_SOURCE_REPLAY','checks':checks,'check_count':len(checks),
 'nonzero_Leibniz_terms':count,'elapsed_seconds':time.monotonic()-start,
 't5_literal_remainder_equality_diagnostic':bool(literal_remainder_equality),
 'scope':'Exact original cubics -> local Taylor matrix -> all t4,t5,t6 coefficients, with only protected denominator/factor divisions. The saturation lift is a separate certificate.'}
save({'checks':checks,'coefficients':coeff},str(out/'pencil_collision_full_source_reader.sobj'))
(out/'pencil_collision_full_source_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'status':report['status'],'checks':len(checks),'seconds':report['elapsed_seconds']},indent=2),flush=True)
