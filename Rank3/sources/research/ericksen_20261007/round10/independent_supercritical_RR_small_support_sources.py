"""Global supercritical R12/H13 sources and literal two-RR comparison.

No D-rank, selected p, det(L), or det(C) inverse is used.  H is the
complete global derivative-minus-source combination, not a derivative
of a previously restricted source.
"""
from sage.all import *
from pathlib import Path
import time,json,os
p=Path(__file__).resolve().parent;start=time.monotonic();checks=[]
def ck(n,v):
    if not bool(v):raise AssertionError(n)
    checks.append(n)
g=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'))
fd=load(str(p/'independent_ratio_weighted_degreefall.sobj'))
bd=load(str(p/'independent_ratio_RR_boundary_leading.sobj'))
hi=load(str(p/'independent_ratio_two_RR_genuine_higher.sobj'))
ex=load(str(p/'independent_ratio_two_RR_mixed_exclusion.sobj'))
ok=g['coefficient_field'];pool=g['polynomial_pool'];K=bd['coefficient_field'];a,r,H,Z=bd['parameters'];q=1+r
G=bd['physical_ring'];T=hi['ring'];v,x,y=T.gens()
ph=ok.ring().hom(bd['parameter_map'],K);cache={}
def cv(c):
    c=ok(c)
    if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
    return cache[c]
def tup(d):return {tuple(map(int,e)):ok(c) for e,c in d.items() if c}
def bi(e):return sum(e[:6]),sum(e[6:])
def selected(d,degree):return {e:c for e,c in tup(d).items() if bi(e)==degree}
def dagface(dag,degree):
    out={}
    for sign,li,ri in dag:
        for e,c in tup(pool[li]).items():
            eb=bi(e)
            for f,d in tup(pool[ri]).items():
                fb=bi(f)
                if (eb[0]+fb[0],eb[1]+fb[1])!=degree:continue
                ef=tuple(e[j]+f[j] for j in range(9));out[ef]=out.get(ef,ok.zero())+sign*c*d
    return {e:c for e,c in out.items() if c}
def mapped(d):return G({e[6:]+e[:6]:cv(c) for e,c in d.items()})
def rv(f,sq):
    dd={}
    for ee,c in T(f).dict().items():
        e=tuple(map(int,ee));ne=(e[0]%2,e[1],e[2]);dd[ne]=dd.get(ne,K.zero())+c*sq**(e[0]//2)
    return T({e:c for e,c in dd.items() if c})
R12=[];H13=[]
for j,it in enumerate(fd['items']):
    raw=tup(it['full'].dict());ck('full R identity '+str(j),raw==tup(pool[g['full_source_refs'][j]]))
    d=selected(raw,(1,2));R12.append(d)
    # The old weighted top is exactly R12+R31.
    old=tup(bd['original_R_top_dictionaries'][j][1]);ck('R top decomposition '+str(j),old=={e:c for e,c in raw.items() if bi(e) in ((1,2),(3,1))})
for j,rec in enumerate(g['records']):
    d=dagface(rec['full_signed_product_DAG'],(1,3));H13.append(d)
    ck('H13 bidegree '+str(j),all(bi(e)==(1,3) for e in d))
    print('global H13',j,'terms',len(d),'seconds',round(time.monotonic()-start,2),flush=True)
items=[]
for m,item in enumerate(hi['items']):
    ch=item['chart'];ev=G.hom(list(ch['RR'])+list(ch['p']),T);sq=ch['v_squared'];values=[]
    for j,d in enumerate(R12):
        value=rv(ev(mapped(d)),sq);ck('R12 vanishes on complete mixed locus %d %d'%(m,j),value==0)
    for j,d in enumerate(H13):
        rec=g['records'][j];i=int(rec['direction']);pos=int(rec['source_position'])
        value=rv(ev(mapped(d)),sq);expected=hi['items'][m]['genuine_R_derivatives'][pos]['values'][i]
        ck('literal global H13 equals genuine weighted derivative %d %d'%(m,j),value==expected)
        values.append(value)
    for ci,cert in enumerate(ex['items'][m]['source_amplitude_certificates']):
        pos=list(g['source_order']).index(int(cert['source_name'][1:]))
        picked=[values[i*6+pos] for i in cert['directions']]
        ck('global H13 amplitude identity %d %d'%(m,ci),sum((f*c for f,c in zip(picked,cert['bezout_coefficients'])),T.fraction_field().zero())==cert['amplitude'])
    items.append({'chart':ch,'H13_values':values,'amplitude_certificates':ex['items'][m]['source_amplitude_certificates']})
payload={'source_files':['independent_ratio_global_R_derivative_degreefall.sobj','independent_ratio_weighted_degreefall.sobj','independent_ratio_RR_boundary_leading.sobj','independent_ratio_two_RR_genuine_higher.sobj','independent_ratio_two_RR_mixed_exclusion.sobj'],'original_field':ok,'coefficient_field':K,'parameter_map':bd['parameter_map'],'original_variable_order':('p12','p13','p21','p23','p31','p32','RR1','RR2','RR3'),'R_order':g['source_order'],'R12_original':R12,'H13_original':H13,'H_order':[(int(it['direction']),int(it['source_position'])) for it in g['records']],'items':items,'checks':checks,'scope':__doc__}
stem='independent_supercritical_RR_small_support_sources';tmp=p/(stem+'.tmp.sobj');save(payload,str(tmp),compress=False);os.replace(tmp,p/(stem+'.sobj'))
(p/(stem+'.json')).write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-start,'R12_terms':[len(d) for d in R12],'H13_terms':[len(d) for d in H13],'scope':__doc__},indent=2)+'\n')
print('PASS',len(checks),'seconds',time.monotonic()-start,flush=True)
