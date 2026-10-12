from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
src=load(str(p/'independent_delta_subcritical_Q0_full_source.sobj'));O=src['ring'];K=O.base_ring();a,r,H,Z=K.gens();al=(Z-H)/(Z*(H+1));be=H*(Z-1)/(Z*(H+1));full=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));old=full['coefficient_field'];pool=full['polynomial_pool'];ph=old.ring().hom([a,a*(1+r)/(1+a*r),H/r,(1+r)*Z/r],K);cache={}
def cv(c):
 if c not in cache:cache[c]=ph(c.numerator())/ph(c.denominator())
 return cache[c]
items=[]
for ci in range(3):
 d=load(str(p/('independent_subcritical_two_RR_full_normal_case'+str(ci)+'.sobj' if ci<2 else 'independent_subcritical_two_RR_full_normal.sobj')));P=d['ring'];t,w,s,e=P.gens();PF=P.fraction_field();vals=list(d['p_images'])+list(d['RR_images']);free=[2,0,1][ci];active=[j for j in range(3) if j!=free]
 flows=[{5:-1/al},{3:-1/(al*t)},{7:1}] if ci==0 else ([{6:1},{1:1/(al*t)},{4:-1/be}] if ci==1 else [{0:be/(al*t)},{6:1},{2:1}])
 def pull(dd,der=None,orig=False):
  out=P.zero()
  for ex,c in dd.items():
   ee=list(map(int,ex));mul=1
   if der is not None:
    mul=ee[der]
    if not mul:continue
    ee[der]-=1
   mon=prod(vals[k]**ee[k] for k in range(9))
   if mon:out+=mul*(cv(c) if orig else K(c))*mon
  return out
 def rawface(f,i,j):return {tuple(map(int,ex)):c for ex,c in f.dict().items() if sum(ex[:6])==i and sum(ex[6:])==j}
 for j,f in enumerate(src['full_R_sources']):
  top=rawface(f,3,1);coeff=[sum(PF(v)*pull(top,k) for k,v in flow.items()) for flow in flows];const=pull(rawface(f,4,0));delta=pull(rawface(f,1,2));ck('R exact critical constant '+str((ci,j)),const==0);ck('R complete critical coefficients '+str((ci,j)),all(coeff[active[k]]==d['matrix'][j,k] for k in range(2)) and coeff[free]==0 and delta==d['matrix'][j,3])
 f=src['full_Q0'];top=rawface(f,1,2);coeff=[sum(PF(v)*pull(top,k) for k,v in flow.items()) for flow in flows];ck('Q0 critical constant',pull(rawface(f,2,1))==0);ck('Q0 complete critical coefficients',all(coeff[active[k]]==d['matrix'][6,k] for k in range(2)) and coeff[free]==0 and d['matrix'][6,3]==0)
 hj=[6,1,0][ci];layers={(3,2):{},(4,1):{},(1,3):{}}
 for sign,ii,jj in full['records'][hj]['full_signed_product_DAG']:
  for ex,c in pool[ii].items():
   for ey,b in pool[jj].items():
    ee=tuple(int(x)+int(y) for x,y in zip(ex,ey));key=(sum(ee[:6]),sum(ee[6:]))
    if key in layers:layers[key][ee]=layers[key].get(ee,old.zero())+sign*c*b
 coeff=[sum(PF(v)*pull(layers[(3,2)],k,True) for k,v in flow.items()) for flow in flows];const=pull(layers[(4,1)],orig=True);delta=pull(layers[(1,3)],orig=True);ck('genuine full-H critical constant '+str(ci),const==0);ck('genuine full-H critical coefficients '+str(ci),all(coeff[active[k]]==d['matrix'][7+hj,k] for k in range(2)) and coeff[free]==0 and delta==d['matrix'][7+hj,3]);items.append({'case':ci,'flows':flows,'H_index':hj,'H_active_coefficients':[coeff[j] for j in active],'H_delta':delta})
# Independently replay all explicit QQ source lifts.
s=load(str(p/'independent_critical_two_RR_R_scaled.sobj'));l=load(str(p/'independent_critical_two_RR_direction_certificate.sobj'))
for it,orig in zip(l['items'],s['items']):
 ck('lift exact source match',it['sources']==orig['QQ_sources'])
 for z in it['items']:ck('full direction source identity',sum(a*b for a,b in zip(z['cofactors'],it['sources']))==z['target'])
pp=load(str(p/'independent_critical_two_RR_direction_H_probe.sobj'));facts=[]
for it in pp['items']:
 u=it['ring'].gen();f=it['first_nonzero'][0];de=it['matrix'].matrix_from_rows(f['rows']).det();ck('actual source three-minor',de==f['determinant']);rem=de%(u*u-it['square']);ck('original degree-two reduction',rem==f['remainder']);ck('det independent of square-root sign',rem.degree()==0);fac=list(K(rem[0]).factor());facts.append(fac)
 for z,exp in fac:
  if z in QQ:continue
  ck('only protected material or pure-base factor',z.degree(K.ring().gen(0))==0 and ((z.degree(K.ring().gen(1))==0) or z/r in QQ or z/(1+r) in QQ))
 ck('nonzero determinant',bool(rem[0]))
 print('CASE',it['case'],'DET FACTOR',str(K(rem[0]).factor()),flush=True)
record={'checks':len(checks),'seconds':time.monotonic()-st,'scope':'All three critical two-RR whole-opposite-pair supports excluded on an explicit material-uniform independent-base open; complete original boundary sources and complex directions retained.'};save({'checks':checks,'items':items,'determinant_factors':facts,'record':record},str(p/'independent_critical_two_RR_source_reader.sobj'),compress=False);(p/'independent_critical_two_RR_source_reader.json').write_text(json.dumps(record,indent=2)+'\n');print('PASS',record,flush=True)
