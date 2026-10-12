"""Exact sign inspection of the only non-protected second-face factor."""
from sage.all import *
from pathlib import Path
import time,json,signal,sys
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(90)
tag=sys.argv[1] if len(sys.argv)>1 else 'plus'
d=load(str(O/f'rr_null_physical_axis_binary_resultant_{tag}.sobj'))
big=max(list(d['factorization']),key=lambda z:z[0].total_degree())[0]
P=PolynomialRing(QQ,names=('s','t'));s,t=P.gens()
re=P({m:QQ(c[0]) for m,c in big.dict().items()});im=P({m:QQ(c[1]) for m,c in big.dict().items()})
V=PolynomialRing(QQ,names=('u','v'));u,v=V.gens();ev=P.hom([1+u,1+u+v],V)
out={};report={}
for label,f in [('real',re),('imag',im)]:
 fac=f.factor();shift=ev(f);co=shift.coefficients()
 report[label]={'factorization':str(fac),'shift_terms':len(co),'positive':sum(z>0 for z in co),'negative':sum(z<0 for z in co)}
 print(label,report[label],flush=True)
 out[label]={'polynomial':f,'factorization':fac,'ordered_positive_shift':shift}
save({'ring':P,'big_factor':big,'parts':out},str(O/f'rr_null_physical_axis_resultant_sign_{tag}.sobj'))
report['seconds']=time.time()-st
(O/f'rr_null_physical_axis_resultant_sign_{tag}.json').write_text(json.dumps(report,indent=2),encoding='utf8')
