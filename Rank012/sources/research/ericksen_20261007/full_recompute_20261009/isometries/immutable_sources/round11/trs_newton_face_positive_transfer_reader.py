"""Independent exact transfer of the pure-rotation positivity obstruction."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent
wall=load(str(OUT/'trs_algebra_newton_face_resultant.sobj'))
cert=load(str(OUT/'tr_infinite_angle_whole_certificate.sobj'))
K=wall['base'];P=K.ring();v,p,j,I1,I2,C,D,h=P.gens();ii=P.base_ring().gen()
source=cert['F'].parent()
def conj_poly(q):
 return source({mon:co.conjugate() for mon,co in q.dict().items()})
hom=source.hom([v,j,I1,I2,C,D],P)
F=hom(conj_poly(cert['F']));G=hom(conj_poly(cert['G']))
top=wall['resultant'].numerator().coefficient({v:2})
ratio=K(top)/(F*G);num=ratio.numerator();den=ratio.denominator();gg=num.gcd(den);ratio=K(num//gg)/K(den//gg)
assert ratio==-j*j/27
report={'status':'PASS','identity':'[v^2] numerator Res = -j^2/27 conjugate(F) conjugate(G)',
        'scope':'All positive simple spectra, arbitrary real KC and KB, nonzero pitch chart.',
        'positivity_source':'tr_infinite_angle_whole_reader.py (507 exact checks)'}
(OUT/'trs_newton_face_positive_transfer_reader.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
