"""All original critical H equations on the constant direction (s,z)=(1,1).

Both Schur columns are retained.  Delta is not assumed constant over H,Z.
"""
from sage.all import *
from pathlib import Path
import time,json,os
p=Path(__file__).resolve().parent;st=time.monotonic()
d=load(str(p/'independent_delta_critical_consistency_minors.sobj'))
K=PolynomialRing(QQ,names=('pa','pr','pH','pZ','pB','pC')).fraction_field()
a,r,H,Z,B,C=K.gens();old=d['field'];ph=old.ring().hom([a,r,H,Z,K.one(),C,B,C],K)
def cv(f):f=old(f);return ph(f.numerator())/ph(f.denominator())
Q=[cv(d['Q_delta_schur']),cv(d['Q_one_schur'])]
items=[]
for it in d['items']:
 hd=cv(it['H_delta_schur']);h0=cv(it['H_one_schur'])
 items.append({'source_row':it['H_row'],'delta':hd,'constant':h0})
 print('POINT_H',it['H_row'],'delta_terms',len(hd.numerator().dict()),'constant_terms',len(h0.numerator().dict()),'seconds',time.monotonic()-st,flush=True)
save({'source_file':'independent_delta_critical_consistency_minors.sobj','field':K,
      'C':d['C'].apply_map(cv),'C_determinant':cv(d['C_determinant']),
      'Q_schur':Q,'items':items,'scope':__doc__},str(p/'independent_delta_critical_point11_schur.sobj'),compress=False)
(p/'independent_delta_critical_point11_schur.json').write_text(json.dumps({'Q_zero':[v==0 for v in Q],'H_constant_zero':[it['constant']==0 for it in items],'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('QZERO',[v==0 for v in Q],flush=True)
