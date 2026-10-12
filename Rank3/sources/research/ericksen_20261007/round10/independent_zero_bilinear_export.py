"""Export three explicit guard identities for CAS-free integer replay."""
from pathlib import Path
import json
from sage.all import load,lcm,ZZ
out=Path(__file__).resolve().parent;d=load(str(out/'independent_zero_bilinear_blocks.sobj'))
def enc(p):
    assert all(c.denominator()==1 for c in p.coefficients())
    return [[list(map(int,mm)),int(c)] for mm,c in sorted(p.dict().items()) if c]
entries=[]
for item in d['data']:
    mult=item['multipliers'];D=lcm([c.denominator() for p in mult for c in p.coefficients()]);assert D>0
    assert item['power']==1
    R=item['ring'];H,Z,S=R.gens()
    gf=[H,Z,H-Z,Z-1,H+1,H-Z+1,S-1,H*S-Z+1,H*S-Z+S]
    removed_decomposition=[]
    for c in item['removed']:
        rest=c;exponents=[]
        for factor in gf:
            e=0
            while rest:
                quotient,remainder=rest.quo_rem(factor)
                if remainder:break
                rest=quotient;e+=1
            exponents.append(e)
        assert rest.is_constant() and rest!=0 and ZZ(rest)==rest
        removed_decomposition.append({'constant':int(rest),'guard_exponents':exponents})
    entries.append({'pair':list(item['pair']),'denominator':int(D),'power':1,
                    'guard':enc(item['guard']),'guard_factors':[enc(p) for p in gf],
                    'sources':[enc(p) for p in item['polynomials']],
                    'multipliers':[enc(D*p) for p in mult],
                    'raw_sources':[enc(p) for p in item['raw']],
                    'removed_factors':[enc(p) for p in item['removed']],
                    'removed_decomposition':removed_decomposition})
data={'variables':['H','Z','S'],'identity':'sum multipliers[i]*sources[i] = denominator*guard',
      'raw_identity':'sources[i]*removed_factors[i] = raw_sources[i]',
      'entries':entries,'scope':'Pointwise zero-block exclusion on the nine-guard chart; no rank-two assertion.'}
(out/'independent_zero_bilinear_integer_certificate.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
print('exported denominators',[entry['denominator'] for entry in entries],flush=True)
