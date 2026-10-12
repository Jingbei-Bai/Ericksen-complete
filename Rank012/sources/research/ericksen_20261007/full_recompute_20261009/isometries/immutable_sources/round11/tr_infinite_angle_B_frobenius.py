from sage.all import *
from pathlib import Path
import time,resource
OUT=Path(__file__).resolve().parent;START=time.time()
f=load(str(OUT/'tr_pure_rotation_logU_frobenius.sobj'))
g=load(str(OUT/'tr_pure_rotation_L_frobenius.sobj'))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
P=f['P'];K=P.fraction_field();v,j,A,T,B,L=map(K,P.gens())
out={}
for label in f['coefficients']:
 z=K(g['coefficients'][label])+(v+L*L)/L*K(f['coefficients'][label])
 n=z.numerator();d=z.denominator();gg=n.gcd(d);z=K(n//gg)/K(d//gg)
 out[label]=z
 print(label,len(z.numerator().dict()),len(z.denominator().dict()),time.time()-START,flush=True)
save({'P':P,'coefficients':out,'note':'L Frobenius +(v+L²)/L logU Frobenius = (v+L²)/(2BL) B Frobenius'},str(OUT/'tr_pure_rotation_Bcompat_frobenius.sobj'))
