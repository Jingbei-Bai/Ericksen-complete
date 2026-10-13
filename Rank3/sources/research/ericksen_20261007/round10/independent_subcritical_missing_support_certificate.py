from sage.all import *
from pathlib import Path
from itertools import product,combinations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
# Combinatorial support exhaustion, without actual zero-pair tangency.
allowed=[]
for mask in product((0,1),repeat=6):
 if not all(mask[2*i] or mask[2*i+1] for i in range(3)):continue
 if all(mask):continue
 if bool(mask[0]*mask[3]*mask[4]) != bool(mask[1]*mask[2]*mask[5]):continue
 allowed.append(mask)
ck('exactly twelve active-D missing-p Delta supports',len(allowed)==12)
torus=[m for m in allowed if m[0]==m[5] and m[1]==m[3] and m[2]==m[4]]
ck('exactly three full-RR supports',len(torus)==3 and all(sum(m)==4 for m in torus))
remaining=[m for m in allowed if m not in torus]
opp=((3,5),(1,4),(0,2))
ck('other nine supports exhausted by three opposite-pair families',all(any(not m[i] and not m[j] for i,j in opp) for m in remaining))
# Six ORIGINAL transitive chart coefficient identities retain z=0.
d=load(str(p/'independent_delta_subcritical_Q0_six_directions.sobj'));P=d['coefficient_ideal_ring'];x,z=P.gens();chart_records=[]
for it in d['items']:
 co=list(it['base_coefficients'].values());cf=it['coefficient_cofactors'][1];target=it['coefficient_targets'][1]
 ck('six-chart explicit original coefficient identity '+str(it['chart']),sum(a*b for a,b in zip(cf,co))==target)
 ck('z-zero exact unit '+str(it['chart']),sum(a.subs({z:0})*b.subs({z:0}) for a,b in zip(cf,co))==-1)
 ck('one nonzero protected material monomial '+str(it['chart']),len(it['material_exponents'])==1)
 chart_records.append({'chart':it['chart'],'coefficient_sources':co,'cofactors':cf,'target':target,'directions':it['directions'],'support_selection':it['support_selection']})
# Every C/Q0 active-e rank-drop condition has an explicit t certificate.
e=load(str(p/'independent_subcritical_two_RR_even_probe.sobj'));R=e['coefficient_ring'];ru,rw,rt=R.gens();ranks=[]
for it in e['items']:
 co=it['base_coefficient_sources'];I=R.ideal(co);cf=rt.lift(I);ck('rankdrop t coefficient identity '+str(it['missing_RR']),sum(a*b for a,b in zip(cf,co))==rt)
 # Rebuild base coefficients from the actual six 2-minors, tracking material powers.
 S=e['direction_ring'];K=e['field'];base=[];matpowers=[]
 for mn in it['minors']:
  pol=S(mn);den=lcm([c.denominator() for c in pol.coefficients()]);cc={};powers=set()
  for ee,c in pol.dict().items():
   val=K(c)*den;ck('minor coefficient denominator clear',val.denominator()==1)
   for ex,b in val.numerator().dict().items():
    ex=tuple(map(int,ex));powers.add(ex[:2]);cc[ex]=cc.get(ex,R.zero())+R({tuple(map(int,ee)):b})
  ck('single material monomial per minor '+str((it['missing_RR'],len(matpowers))),len(powers)<=1)
  matpowers.append(powers);base+=list(cc.values())
 ck('all original minor coefficient sources exact '+str(it['missing_RR']),base==co)
 ranks.append({'missing_RR':it['missing_RR'],'source_coefficients':co,'target':rt,'cofactors':cf,'material_powers':matpowers})
save({'checks':checks,'allowed_supports':allowed,'torus_supports':torus,'two_RR_supports':remaining,'six_chart_unit_certificates':chart_records,'rankdrop_certificates':ranks},str(p/'independent_subcritical_missing_support_certificate.sobj'),compress=False)
(p/'independent_subcritical_missing_support_certificate.json').write_text(json.dumps({'checks':len(checks),'seconds':time.monotonic()-st,'rankdrop_cofactor_nonzero':[sum(bool(f) for f in it['cofactors']) for it in ranks]},indent=2)+'\n');print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
