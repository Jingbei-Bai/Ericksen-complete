"""Bind the complete partial gauge to the manuscript state derivation."""
from sage.all import *
from pathlib import Path
import hashlib,json
HERE=Path(__file__).absolute().parent
INPUT=HERE/'source_inputs'
d=load(str(INPUT/'independent_ratio_p_boundary_second_gauge.sobj'))
g=load(str(INPUT/'independent_partial_D_Newton_gauges_reader.sobj'))
R=g['ring'];K=R.base_ring();v=R.gens()
cp=lambda f:R({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
hat=matrix(R,3,9,[cp(f) for f in d['first_gauged_fields'].list()])
base=matrix(R,3,2,[cp(f) for f in d['base_fields'].list()])
src=[cp(f) for f in d['six_full_sources']]
item=g['items'][2]
assert item['direction']==1 and item['source_row']==2 and item['p_index']==2 and item['source_sign']==1
expected=vector(R,[v[2]*hat[1,j] for j in range(9)])
unit=K(d['source_unit_coefficients'][1])
expected[7]-=unit*src[2]
assert expected==item['field_numerators']
assert vector(R,[v[2]*base[1,j] for j in range(2)])==item['base_numerators']
J=FractionField(PolynomialRing(QQ,names=('k','r','eta','zeta')))
k,r,eta,zeta=J.gens();q=1+r;m=1+k
ph=K.ring().hom([k/r,k*q/(r*m),eta/r,q*zeta/r],J)
coefficient=ph(unit.numerator())/ph(unit.denominator())
assert coefficient==-2*(zeta-1)/zeta,coefficient

report={'status':'MANUSCRIPT_OPERATOR_MATCH',
        'operator':'p21*Ehat2-c2diag*R2*d/dR2',
        'c2diag':'-2*(zeta-1)/zeta','ratio_components_match':True,
        'input_hashes':{n:hashlib.sha256((INPUT/n).read_bytes()).hexdigest()
         for n in ['independent_ratio_p_boundary_second_gauge.sobj',
                   'independent_partial_D_Newton_gauges_reader.sobj']}}
(HERE/'manuscript_gauge_binding.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
