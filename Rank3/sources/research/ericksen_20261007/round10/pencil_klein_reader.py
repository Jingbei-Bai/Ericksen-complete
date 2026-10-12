"""Exact parity compression only; no coprimality claim or dropped wall."""
from pathlib import Path
import json,time
from sage.all import load,save,PolynomialRing,QQ
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_fast_numerators.sobj'))
P=PolynomialRing(QQ,names=('U','V','S','a','q','T','Hb'))
A={};B={};C={};D={}
for m,co in d['rows'][0].dict().items():
    par=tuple(e%2 for e in m[:3]);assert par in [(0,0,1),(1,1,0)]
    target=A if par==(0,0,1) else B
    mon=tuple((m[i]-par[i])//2 for i in range(3))+tuple(m[3:]);target[mon]=co
for m,co in d['rows'][1].dict().items():
    par=tuple(e%2 for e in m[:3]);assert par in [(0,0,0),(1,1,1)]
    target=C if par==(0,0,0) else D
    mon=tuple((m[i]-par[i])//2 for i in range(3))+tuple(m[3:]);target[mon]=co
arr=list(map(P,(A,B,C,D)))
F=d['flat_ring']
def reconstruct(parts):
    result={}
    for poly,shift in parts:
        for m,co in poly.dict().items():
            mon=tuple(2*m[i]+shift[i] for i in range(3))+tuple(m[3:])
            assert mon not in result
            result[mon]=co
    return F(result)
checks={'H_parity_identity':d['rows'][0]==reconstruct([(arr[0],(0,0,1)),(arr[1],(1,1,0))]),
        'Q0_parity_identity':d['rows'][1]==reconstruct([(arr[2],(0,0,0)),(arr[3],(1,1,1))])}
assert all(checks.values())
records=[{'name':name,'invariant_degree':max(sum(m[:3]) for m in f.dict()),
          'invariant_support':len(set(tuple(m[:3]) for m in f.dict())),
          'flat_terms':len(f.dict())} for name,f in zip('ABCD',arr)]
save({'ring':P,'rows':arr,'checks':checks,'records':records,
      'convention':'H=s*A(U,V,S)+u*v*B(U,V,S), Q0=C(U,V,S)+u*v*s*D(U,V,S); U=u^2,V=v^2,S=s^2.',
      'scope':'No common factor was removed, and no determinant norm was declared protected.'},str(out/'pencil_klein_reader.sobj'))
(out/'pencil_klein_reader.json').write_text(json.dumps({'status':'TWO_EXACT_PARITY_REPLAYS_PASS','checks':checks,'records':records,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
print(json.dumps({'checks':checks,'records':records,'seconds':time.monotonic()-start},indent=2),flush=True)
