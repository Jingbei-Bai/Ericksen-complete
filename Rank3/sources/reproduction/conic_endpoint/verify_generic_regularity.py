"""Full generic regularity: exact negative-face cancellation, no spectrum specialization."""
from sage.all import *
from pathlib import Path
import time,json
out=Path(__file__).absolute().parent;p=out/'source_inputs'
start=time.monotonic();checks=[]
def ck(label,b):
    assert bool(b),label
    checks.append(label)
qf=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'))
ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'))
R=qf['ring'];ell=qf['L']*vector(R,R.gens()[6:])
def layer(f,di,dj):return R({tuple(map(int,e)):c for e,c in f.dict().items() if sum(e[:6])==di and sum(e[6:])==dj})
DC=matrix(R,qf['D'])*matrix(R,qf['C'])
for i,f in enumerate(qf['full_R_sources']):
    ck('R'+str(i)+' exact 31 L factorization',layer(f,3,1)==sum(DC[i,j]*ell[j] for j in range(3)))
    for e in f.dict():
        di,dj=sum(e[:6]),sum(e[6:])
        if (di,dj)!=(3,1):ck('R regular '+str((i,di,dj)),8-2*di-3*dj>=0 and 4-di-dj>=0)
f=qf['full_Q0']
ck('Q0 exact 12 L factorization',layer(f,1,2)==sum(qf['Q0_L_cofactors'][j]*ell[j] for j in range(3)))
for e in f.dict():
    di,dj=sum(e[:6]),sum(e[6:])
    if (di,dj)!=(1,2):ck('Q0 regular '+str((di,dj)),7-2*di-3*dj>=0 and 3-di-dj>=0)
for item in ev['items']:
    remainder=item['full_source']-item['square_unit']*ell[int(item['square_component'])]**2-sum(item['A_p_quadratic'][j]*ell[j] for j in range(3))
    for e in remainder.dict():
        di,dj=sum(e[:6]),sum(e[6:])
        ck('even regular '+str((item['Q_index'],di,dj)),8-2*di-3*dj>=0 and 4-di-dj>=0)
print('R Q0 even generic regularity PASS',len(checks),flush=True)
hs=load(str(p/'independent_delta_subcritical_H_initials.sobj'))
full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
K=full['coefficient_field'];a,b,h,z=K.gens()
r=(a-b)/(a*(b-1));q=1+r;H=r*h;Z=r*z/q
alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1))
O=PolynomialRing(K,names=full['variable_names']);v=O.gens()
L=[v[0]*v[6]-alpha*v[5]*v[8],v[1]*v[6]-alpha/beta*v[3]*v[7],v[2]*v[7]-beta*v[4]*v[8]]
pool=full['polynomial_pool'];cache={}
def layers(i):
    if i not in cache:
        dd={}
        for e,c in pool[i].items():
            ee=tuple(map(int,e));degree=(sum(ee[:6]),sum(ee[6:]))
            dd.setdefault(degree,{})[ee]=c
        cache[i]={degree:O(f) for degree,f in dd.items()}
    return cache[i]
for i,record in enumerate(full['records']):
    bad={};possible=set()
    for sign,li,ri in record['full_signed_product_DAG']:
        for (di,dj),f in layers(li).items():
            for (ei,ej),g in layers(ri).items():
                degree=(di+ei,dj+ej);possible.add(degree)
                if degree==(3,2) or 11-2*degree[0]-3*degree[1]<0 or 5-sum(degree)<0:
                    bad[degree]=bad.get(degree,O.zero())+sign*f*g
    expected=sum(O({tuple(map(int,e)):c for e,c in hs['items'][i]['L_cofactors_original'][j].items()})*L[j] for j in range(3))
    ck('H'+str(i)+' exact original 32 L factorization',bad.get((3,2),O.zero())==expected)
    for degree,f in bad.items():
        if degree!=(3,2):ck('H'+str(i)+' zero forbidden full layer '+str(degree),f==0)
    print('H',i,'generic L and forbidden-layer cancellation PASS',round(time.monotonic()-start,2),flush=True)
payload={'status':'PASS','scope':__doc__,'checks':len(checks),'identities':checks,'seconds':time.monotonic()-start}
(out/'generic_regularity.json').write_text(json.dumps(payload,indent=2)+'\n')
print('PASS',len(checks),'seconds',time.monotonic()-start,flush=True)
