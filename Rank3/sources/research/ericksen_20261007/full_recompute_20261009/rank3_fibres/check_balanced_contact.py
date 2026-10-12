from sage.all import *
from pathlib import Path
import json
out=Path(__file__).absolute().parent;d=load(str(out.parent/'foundations/state17_fresh.sobj'))
P=d['ring'];xx=P.gens();a,b=P.base_ring().gens();A=a*(b-1);C=b*(a-1);D=a-b;r=D/A;q=C/A
H=xx[4];Z=C*(A+D*H)/(A*D)
images=list(xx);images[0]=P.one();images[8]=P(Z);ph=P.hom(images,P)
B=[Z*xx[5]+C/A*H*xx[7],xx[6]+D/C*Z*xx[2],H*xx[1]-A/D*xx[3]]
expected=[-2*D/C*B[0],2*B[1],2*D/A*B[2]]
checks={}
for i in range(3):
 f=d['F'][i]
 derivative=ph(r*f[4]-r/q*f[8]+(-r*H+r/q*Z)*f[0])
 checks['balanced_line_actual_contact_'+str(i)]=derivative==expected[i]
assert all(checks.values()),checks
(out/'balanced_contact.json').write_text(json.dumps(dict(status='PASS',checks=checks),indent=2));print(checks)
