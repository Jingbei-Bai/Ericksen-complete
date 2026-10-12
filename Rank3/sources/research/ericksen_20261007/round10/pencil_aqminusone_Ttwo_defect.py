"""Exact local root-lifting defect at the only common initial root."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'pencil_aqminusone_Ttwo.sobj'))
C=PolynomialRing(QQ,'q');q=C.gen();K=C.fraction_field();R=PolynomialRing(K,'X');X=R.gen();series=[]
for p,ht in zip(d['cubics'],(14,17)):
    series.append([R([C({int(mon[0]):v for mon,v in co.dict().items() if mon[1]==ht+j-k}) for j,co in enumerate(p.list())]) for k in range(3)])
F0,F1,F2=series[0];G0,G1,G2=series[1];root=K(-2*(q*q-10*q+1)/(q-1)**3)
assert not F0(root) and not G0(root) and not G0.derivative()(root)
g1=G1(root);fp=F0.derivative()(root);shift=-F1(root)/fp
print('Fprime',fp.factor(),flush=True);print('G1',g1.factor() if g1 else '0',flush=True)
second=G2(root)+shift*G1.derivative()(root)+shift**2*G0.derivative(2)(root)/2
print('SECOND',second.factor() if second else '0',flush=True)
save({'series':series,'root':root,'Fprime':fp,'G1_at_root':g1,'shift':shift,'effective_second':second},str(out/'pencil_aqminusone_Ttwo_defect.sobj'))
(out/'pencil_aqminusone_Ttwo_defect.json').write_text(json.dumps({'Fprime':str(fp.factor()),'G1':str(g1.factor()) if g1 else '0','effective_second':str(second.factor()) if second else '0','elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
