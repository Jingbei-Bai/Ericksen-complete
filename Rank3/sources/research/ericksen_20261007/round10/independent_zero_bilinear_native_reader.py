"""Exact integer proof replay using only the Python standard library."""
from pathlib import Path
import json,time
out=Path(__file__).resolve().parent;start=time.monotonic()
d=json.loads((out/'independent_zero_bilinear_integer_certificate.json').read_text())
def dec(p):return {tuple(mm):int(c) for mm,c in p if c}
def add(p,q):
    p=dict(p)
    for mm,c in q.items():
        v=p.get(mm,0)+c
        if v:p[mm]=v
        else:p.pop(mm,None)
    return p
def mul(p,q):
    ans={}
    for i,c in p.items():
        for j,e in q.items():
            mm=tuple(a+b for a,b in zip(i,j));ans[mm]=ans.get(mm,0)+c*e
    return {mm:c for mm,c in ans.items() if c}
reports=[]
for item in d['entries']:
    D=item['denominator'];assert isinstance(D,int) and D>0 and item['power']==1
    checks={};guard=dec(item['guard']);prod={(0,0,0):1}
    for f in item['guard_factors']:prod=mul(prod,dec(f))
    checks['nine_guard_product']=prod==guard
    sources=[dec(p) for p in item['sources']]
    assert all(sources)
    total={}
    for i,(p,m,raw,c) in enumerate(zip(sources,item['multipliers'],item['raw_sources'],item['removed_factors'])):
        checks['raw_source_'+str(i)]=mul(p,dec(c))==dec(raw)
        decomposition=item['removed_decomposition'][i]
        constant=decomposition['constant'];assert isinstance(constant,int) and constant!=0
        factor={(0,0,0):constant}
        for f,exponent in zip(item['guard_factors'],decomposition['guard_exponents']):
            assert isinstance(exponent,int) and exponent>=0
            for _ in range(exponent):factor=mul(factor,dec(f))
        checks['removed_factor_protected_'+str(i)]=factor==dec(c)
        total=add(total,mul(p,dec(m)))
    checks['explicit_integer_identity']=total=={mm:D*c for mm,c in guard.items()}
    assert all(checks.values()),(item['pair'],checks)
    reports.append({'pair':item['pair'],'checks':checks,'denominator':D})
report={'status':'EXACT_INTEGER_CERTIFICATES_PASS','reports':reports,'elapsed_seconds':time.monotonic()-start}
(out/'independent_zero_bilinear_native_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
