"""Printed graph and conic identities on the cyclic single-range boundary."""
from pathlib import Path
import json,hashlib,time
from sage.all import *
h=Path(__file__).absolute().parent
replay=h/'replay'/'round10'
checks=[];start=time.monotonic()
d=load(str(replay/'independent_critical_single_RR3_normal.sobj'));P=d['ring'];a,r,H,Z=P.base_ring().gens();x,z,w,e0,e1,e2,de=P.gens();s=x+z;q=1+r;k=a*r;m=1+k
f0=-2*x*w+(4*H*Z-2*H-4*Z**2+2*Z)/(H*Z)*s*w
f1=-2*Z*(H-1)/(H*(Z-1))*w
f2=2*x+(2*H*Z-2*H-4*Z**2+2*Z)/Z*s
dd=4*r*(H+1)**3/(q*H**2)*s*w
assert list(d['graph_images'][3:])==[f0,f1,f2,dd]
checks.append('printed appendix four whole graph functions')
C=d['graph_boundary'][27]
printed=-4*k/m**2*x**4+8*k/m**2*x**3*s+8*(k+w*w)/m*x*x-8*(3*k-w*w)/m*x*s
# Polynomial substitution z=s-x and division by the literal square of s.
T=PolynomialRing(P.base_ring(),names=('x','s','w'));xx,ss,ww=T.gens();ph=P.hom([xx,ss-xx,ww,0,0,0,0],T)
delta=ph(C-printed)
assert all(e[1]>=2 for e in delta.dict())
checks.append('printed appendix exact conic C modulo s²')
source=replay/'independent_critical_single_RR3_normal.sobj'
report=dict(status='PASS',checks=checks,source=str(source),
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            seconds=time.monotonic()-start)
(h/'source_bridges.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS',len(checks),'printed boundary identities',flush=True)
