"""First two asymptotic coefficients at the common Newton root."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_vzero_sone_cubics.sobj'));old=d['base_ring']
A=PolynomialRing(QQ,names=('a','q','T'));a,q,T=A.gens();K=A.fraction_field();R=PolynomialRing(K,'X');X=R.gen()
series=[]
for p,height in zip(d['cubics'],(20,23)):
    coeffs=[]
    for k in range(3):
        ans=R.zero()
        for j,co in enumerate(p.list()):
            ex=height-j-k
            val=A({tuple(mon[:3]):cc for mon,cc in co.dict().items() if mon[3]==ex})
            ans+=val*X**j
        coeffs.append(ans)
    series.append(coeffs)
F0,F1,F2=series[0];G0,G1,G2=series[1];root=K(T*(T-2)*(a*q+1)/4)
assert not F0(root) and not G0(root) and not G0.derivative()(root)
g1=G1(root);f1=F1(root)
report={'G1_at_common_root_zero':not bool(g1),'G1_at_common_root_factorization':str(g1.factor()) if g1 else '0',
        'F1_at_common_root_factorization':str(f1.factor()) if f1 else '0'}
print(json.dumps(report),flush=True)
if not g1:
    fprime=F0.derivative()(root);assert fprime
    shift=-f1/fprime
    second=G2(root)+shift*G1.derivative()(root)+shift*shift*G0.derivative(2)(root)/2
    report['G_effective_second_factorization']=str(second.factor()) if second else '0'
    report['G_effective_second_zero']=not bool(second)
    report['common_root_shift_factorization']=str(shift.factor()) if shift else '0'
    data={'series':series,'root':root,'shift':shift,'effective_second':second}
else:data={'series':series,'root':root,'G1_at_root':g1}
save(data,str(out/'pencil_newton_common_root.sobj'))
report['elapsed_seconds']=time.monotonic()-start
(out/'pencil_newton_common_root.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
