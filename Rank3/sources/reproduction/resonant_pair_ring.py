"""Define the polynomial ring used by the complete three-pair pullback."""
from pathlib import Path
from sage.all import QQ,PolynomialRing,save
root=Path(__file__).absolute().parent.parent
K=PolynomialRing(QQ,names=('lf_a','lf_r')).fraction_field()
Q=PolynomialRing(K,names=('lf_l','lf_rho','lf_k','lf_s','lf_w','lf_e0','lf_e1','lf_e2','lf_delta'))
path=root/'reproduction/resonant/independent_inactive_RRtorus_lifted_fields_cyclic.sobj'
path.parent.mkdir(parents=True,exist_ok=True)
save({'ring':Q},str(path),compress=False)
assert Q.ngens()==9 and K.ngens()==2
print('PASS exact resonant pair polynomial ring',flush=True)
