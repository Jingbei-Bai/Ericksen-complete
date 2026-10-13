"""Critical delta-nonzero line consequences from the genuine affine sources."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_delta_critical_consistency_minors.sobj'));old=d['field'];K=PolynomialRing(QQ,names=('la','lr','lt')).fraction_field();a,r,t=K.gens();P=PolynomialRing(K,names=('lH','lZ','lB','lC'));H,Z,B,C=P.gens();ph=old.ring().hom([P(a),P(r),H,Z,P(-t/(t+1)),C*t,B,C],P)
sources=[('Q0',d['Q_delta_schur'])]+[('H'+str(it['H_row']),it['H_delta_schur']) for it in d['items']];rows=[];labels=[];items=[];monsall=set()
for label,f in sources:
 n=ph(f.numerator());den=ph(f.denominator());mins=tuple(min(int(e[k]) for e in n.dict()) for k in (2,3));unit=B**mins[0]*C**mins[1];nn=n//unit;groups={}
 for ee,c in nn.dict().items():
  ex=tuple(map(int,ee));key=ex[:2];groups.setdefault(key,{})[ex[2:]]=c;monsall.add(ex[2:])
 items.append({'label':label,'source':f,'full_numerator_pull':n,'full_denominator_pull':den,'removed_amplitude_unit':unit,'coefficient_dicts':groups});print('LINE_SOURCE',label,'amplitude_support',sorted(set(e[2:] for e in nn.dict())),'seconds',time.monotonic()-st,flush=True)
mons=sorted(monsall)
for it in items:
 for key,dd in sorted(it['coefficient_dicts'].items()):rows.append([dd.get(e,K.zero()) for e in mons]);labels.append((it['label'],key))
M=matrix(K,rows);print('MATRIX',M.nrows(),M.ncols(),'mons',mons,'rank',M.rank(),'seconds',time.monotonic()-st,flush=True)
save({'field':K,'ring':P,'sources':items,'matrix':M,'monomials':mons,'labels':labels,'scope':__doc__},str(p/'independent_delta_critical_line_delta_coefficients.sobj'),compress=False)
(p/'independent_delta_critical_line_delta_coefficients.json').write_text(json.dumps({'rows':M.nrows(),'columns':M.ncols(),'monomials':mons,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
