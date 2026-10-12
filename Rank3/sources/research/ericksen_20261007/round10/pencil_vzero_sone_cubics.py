"""Exact degree-preserving candidate state slice; all base parameters free.

No resultant or coprimality conclusion is asserted by this export.
"""
from pathlib import Path
import time,json
from sage.all import load,save,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_state_slice_probe.sobj'));CR=d['base_ring'];K=CR.fraction_field();a,q,T,Hb=map(K,CR.gens())
rows=d['coefficient_dicts_by_u_power'];S=PolynomialRing(CR,names=('U',));U=S.gen();slices=[];records=[]
for i,row in enumerate(rows):
 dic={}
 for up,terms in row.items():
  c=sum((co for (vp,sp),co in terms.items() if vp==0),K.zero())
  if c:
   assert up%2==0
   assert c.denominator()==1
   dic[(up//2,)]=CR(c.numerator())
 f=S(dic);slices.append(f)
 rec={'source':i,'degree_U':int(f.degree()),'coefficient_base_terms':[c.number_of_terms() for c in f.list()]}
 records.append(rec);print(rec,round(time.monotonic()-start,2),flush=True)
c6=rows[0][6][(0,1)];d62=rows[1][6][(0,2)];d60=rows[1][6][(0,0)]
checks={'H_u6_is_c6_s':set(rows[0][6])=={(0,1)},
 'Q0_u6_is_d62_s2_plus_d60':set(rows[1][6])=={(0,2),(0,0)},
 'H_slice_degree6':slices[0].degree()==3,'Q0_slice_degree6':slices[1].degree()==3,
 'H_slice_leading':slices[0].leading_coefficient()==c6,
 'Q0_slice_leading':slices[1].leading_coefficient()==d62+d60}
assert all(checks.values()),checks
oldguard=CR(a*q*(q-1)*(1+a*(q-1))*T*(T-1)*Hb*(Hb+T)*(Hb+1)*(1-Hb-T)*(a-1))
payload={'base_ring':CR,'slice_ring':S,'parameters':tuple(map(K,CR.gens())),
 'cubics':slices,'c6':c6,'d62':d62,'d60':d60,'q0_slice_leading':d62+d60,
 'old_guard':oldguard,'checks':checks,'records':records,
 'slice':'v=0,s=1,U=u^2; source 0=primitive H, source 1=primitive Q0',
 'scope':'Exact state slice. Uniform generic-Hb nonzero leading coefficients and resultant must still be proved separately.'}
save(payload,str(out/'pencil_vzero_sone_cubics.sobj'))
print('saved_slice',round(time.monotonic()-start,2),flush=True)
factor_records={}
for name,co in [('c6',c6),('q0_leading',d62+d60)]:
 facs=[]
 for f,e in CR(co).factor():
  facs.append({'factor':str(f),'power':int(e),'old_guard':not oldguard.quo_rem(f)[1],
               'degree_Hb':int(f.degree(CR.gen(3)))})
 factor_records[name]=facs
 print(name,[(x['power'],x['old_guard'],x['degree_Hb'],len(x['factor'])) for x in facs],round(time.monotonic()-start,2),flush=True)
payload['leading_factor_records']=factor_records
save(payload,str(out/'pencil_vzero_sone_cubics.sobj'))
report={'status':'PASS_EXACT_ALL_PARAMETER_STATE_SLICE','checks':checks,'records':records,
 'leading_factor_records':factor_records,'elapsed_seconds':time.monotonic()-start,'scope':payload['scope']}
(out/'pencil_vzero_sone_cubics.json').write_text(json.dumps(report,indent=2)+'\n')
print('DONE',round(time.monotonic()-start,2),flush=True)
