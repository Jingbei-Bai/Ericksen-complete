"""Independent multiplyback of coupled minors and the all-material V lift.

Keeps every row content, the eta/qk pivots and material-amplitude walls.
The geometric quotient/differential argument is recorded separately.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(s,b):assert b,s;checks.append(s)
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+900*1024**2,)*2)
d=load(str(p/'independent_ordinary_minus_two_coupled_minor.sobj'));hc=load(str(p/'independent_ordinary_minus_two_horizontal_coefficients.sobj'));cert=load(str(p/'independent_ordinary_minus_two_V_certificate.sobj'));orig=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));quad=load(str(p/'independent_incoming_simple_wall_quadratics.sobj'))
P=d['ring'];a,r,H,Z,l,u=P.gens();U=orig['ring'];ol,ou,ov=U.gens();K=d['coefficient_map'].codomain();ma,mr,mh,mz=K.gens();ph=d['coefficient_map'];Kp=K.ring();M=orig['matrix']
def toP(f,den):
 out={}
 for ee,co in U(f).dict().items():
  e=tuple(map(int,ee));ck('no remaining v exponent '+str(len(checks)),e[2]==0)
  for be,bc in Kp(den*ph(co)).dict().items():
   key=tuple(map(int,be))+(e[0],e[1]);out[key]=out.get(key,0)+bc
 return P(out)
for j,row in d['cleared_rows'].items():
 got=[U(M[j,c].subs({ov:-2})) for c in (0,1,2,5)]+[U((M[j,3]-M[j,0]).derivative(ov).subs({ov:-2}))]
 ck('whole v=-2 kernel row '+str(j),(M[j,3]-M[j,0]).subs({ov:-2})==0);ck('original five columns row '+str(j),got==d['original_rows'][j]);den,ct=d['row_scalings'][j]
 for k,f in enumerate(got):ck('exact content+den multiplyback %d %d'%(j,k),toP(f,den)==ct*row[k])
 ck('rowcontent retains only state l '+str(j),ct.degree(u)==0)
for it in d['items']:
 rows=it['rows'];A=matrix(P,[d['cleared_rows'][j] for j in rows]);ck('primitive determinant exact '+str(rows),A.det(algorithm='df')==it['primitive_determinant']);ck('factorization exact '+str(rows),P(it['factorization'])==it['primitive_determinant'])
 others=[f for f,e in it['factorization'] if not(f.degree(l)>1 or f.degree(u)>0)]
 print('OUTER FACTORS',rows,[(str(f),e) for f,e in it['factorization'] if f in others],flush=True)
Q=hc['ring'];C=Q.base_ring();ca,cr,cz=C.gens();L,V=Q.gens();allc=[]
for it,ref in zip(d['items'],hc['items']):
 state=[f for f,e in it['factorization'] if (f.degree(l)>1 or f.degree(u)>0) and len(f.dict())>1];ck('one retained statefactor '+str(it['rows']),len(state)==1);f=state[0];byH={}
 for ee,co in f.dict().items():
  e=tuple(map(int,ee));ck('even l/u powers '+str(len(checks)),e[4]%2==e[5]%2==0);hh=e[2]+e[5]//2;val=C(co)*ca**e[0]*cr**e[1]*cz**e[3]*L**(e[4]//2)*V**(e[5]//2);byH[hh]=byH.get(hh,Q(0))+val
 ck('exact H coefficients '+str(it['rows']),byH==ref['coefficients']);ck('degree4 H cover '+str(it['rows']),set(byH)==set(range(5)));allc.extend(byH[h] for h in sorted(byH))
ck('all fifteen original coefficients',allc==hc['all_coefficients']==cert['sources']);cofs=cert['cofactors'];ck('exact V cofactor multiplyback',sum((f*c for f,c in zip(allc,cofs)),Q(0))==V);ck('cofactors state constant',all(not f or f.total_degree()==0 for f in cofs))
den=lcm([C(co).denominator() for f in cofs for co in f.dict().values()]);ck('all denominator factors retained',den==cert['denominator']);zp=C.ring().gen(2)
for it in cert['denominator_factors']:
 f=it['factor'];dg=f.degree(zp);lead=f.coefficient({zp:dg});ck('den base degree+leading exact '+str(len(checks)),dg==it['Z_degree'] and lead==it['Z_leading'])
 if dg==12:ck('all-material protected large factor',lead==66816*(ca*cr+1)**2*(cr+1)**2)
 elif dg==0:ck('only a/r pure material factors',f in (C.ring().gen(0),C.ring().gen(1)))
 else:ck('monic pure base factors',lead==1)
ck('eta pivot comes from fullsource',d['eta_pivot']==M[0,4]);qk=next(x for x in quad['items'] if x['name']=='k');ck('qk2 pivot comes from fullsource',d['quadratic_pivot']==qk['values'][2]['second'])
for it in d['items']:
 selected=[0,2]+list(it['rows'])
 for j in selected:
  if j!=0:ck('full seven-column eta sparse row '+str(j),M[j,4].subs({ov:-2})==0)
  if j!=2:ck('full seven-column raw quadratic sparse row '+str(j),qk['values'][j]['second']==0)
 ck('eta pivot nonzero symbolic '+str(it['rows']),M[0,4]!=0)
 ck('quadratic pivot nonzero symbolic '+str(it['rows']),qk['values'][2]['second']!=0)
 print('FULL7 LAPLACE VERIFIED',it['rows'],flush=True)
 etaunit=ph(d['eta_pivot'].monomial_coefficient(ou));print('ETA COEFFICIENT TRANSFORMED',etaunit.factor(),flush=True)
 qp=U(d['quadratic_pivot']);qmapped={}
 for ee,co in qp.dict().items():qmapped[tuple(map(int,ee))]=ph(co)
 print('QUADRATIC PIVOT TRANSFORMED',qmapped,flush=True)
rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'all_material_denominator_unit':True,'retained_walls':['l=0','u=0','m*l^2=k'],'scope':__doc__}
save({'record':rec,'checks':checks,'denominator':den,'cofactors':cofs,'sources':allc},str(p/'independent_ordinary_minus_two_V_reader.sobj'),compress=False);(p/'independent_ordinary_minus_two_V_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,flush=True)
