"""Reconstruct ALL 22 critical affine sources and all three even quadrics.

Source coefficients come from fresh Section2 joint fields and complete H DAG.
Polynomial module cofactors are explicit imported witnesses, each multiplied
against the fresh source before use. No zero remainder report is trusted.
"""
from sage.all import *
from pathlib import Path
import time,json,sys,hashlib
out=Path(__file__).absolute().parent;research=out.parent.parent;st=time.monotonic();checks=[]
partial=bool(len(sys.argv)>1);bounds=list(map(int,sys.argv[1:3])) if partial else [0,18]
rowout=out/'critical_rows';rowout.mkdir(exist_ok=True)
f=load(str(out/'fresh_global_polynomials.sobj'));joint=load(str(out/'fresh_joint_fields.sobj'));dag=load(str(out/'fresh_global_source_DAG.sobj'))
OR=f['ring'];OK=OR.base_ring();pool=[OR(dd) for dd in dag['polynomial_pool']]
K=PolynomialRing(QQ,names=('a','r','H','Z','p21','p23','p31','p32')).fraction_field();a,r,H,Z,s,t,B,C=K.gens();q=1+r;m=1+a*r
P=PolynomialRing(K,names=('e0','e1','e2'));ee=P.gens();alpha=(Z-H)/(Z*(H+1));beta=H*(Z-1)/(Z*(H+1));vals=(K.one(),t*B/(s*C),s,t,B,C,alpha*C,beta*B/s,K.one())
phi=OK.ring().hom([a,a*q/m,H/r,q*Z/r],K);cache={}
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=phi(c.numerator())/phi(c.denominator())
 return cache[c]
def pull(poly):return sum((P(cv(c))*prod(vals[j]**int(e[j]) for j in range(9)) for e,c in poly.dict().items()),P.zero())
def part(poly,bi):return OR({tuple(e):c for e,c in poly.dict().items() if (sum(e[:6]),sum(e[6:]))==bi})
def row(linear,delta):
 ck('affine normal source '+str(len(checks)),linear.total_degree()<=1 and delta.total_degree()<=0)
 return [K(linear.derivative(e)) for e in ee]+[K(delta),K(linear.constant_coefficient())]
C0=matrix(P,[[pull(v) for v in row] for row in f['C'].rows()]);o=vector(P,[pull(v) for v in f['offset']]);R0=C0*vector(P,ee)+o;rows=[]
for j,(idx,unit) in enumerate([(0,B),(2,s),(3,K.one())]):rows.append(row(R0[j],pull(part(f['sources'][idx],(1,2)))/unit))
Q12=part(f['Q0'],(1,2));qb=[]
for pj,rj in [(0,0),(1,0),(2,1)]:
 dd={}
 for e,c in Q12.dict().items():
  e=list(e)
  if e[pj]:e[pj]-=1;e[6+rj]-=1;ck('Q0 cofactor exponents',min(e)>=0);dd[tuple(e)]=c
 qb.append(OR(dd))
ell=f['L']*vector(OR,OR.gens()[6:]);ck('fresh exact Q0 module',sum(qb[j]*ell[j] for j in range(3))==Q12)
rows.append(row(sum(pull(qb[j])*ee[j] for j in range(3))+pull(part(f['Q0'],(2,1))),P.zero()))
# Witnesses contain the normalized polynomial coefficients, never just a flag.
hm=load(str(research/'round10/independent_delta_H_global_L_module.sobj'));scales=[OK(v) for v in hm['original_to_binomial_variable_scales']]
fingerprints={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [out/'fresh_global_polynomials.sobj',out/'fresh_global_source_DAG.sobj',research/'round10/independent_delta_H_global_L_module.sobj']}
hmK=hm['original_field'];hphi=hmK.ring().hom(list(OK.gens()),OK)
def hc(c):c=hmK(c);return hphi(c.numerator())/hphi(c.denominator())
def cofactor(poly):return OR({tuple(map(int,e)):hc(c)/prod(scales[k]**int(e[k]) for k in range(9)) for e,c in poly.dict().items()})
savedH=[]
for n,rec in enumerate(dag['records']):
 if partial and not bounds[0]<=n<bounds[1]:continue
 checkpoint=rowout/('row'+str(n)+'.sobj')
 if checkpoint.is_file():
  cached=load(str(checkpoint));ck('fresh row input hashes '+str(n),sorted(cached['fingerprints'].values())==sorted(fingerprints.values()));checks.extend(cached['checks']);rows.append(cached['coefficients']);savedH.append(cached);continue
 begin_checks=len(checks)
 pieces={(3,2):{},(4,1):{},(1,3):{}}
 for sign,li,ri in rec['full_signed_product_DAG']:
  for ex,co in pool[li].dict().items():
   for ey,co2 in pool[ri].dict().items():
    e=tuple(int(x)+int(y) for x,y in zip(ex,ey));bi=(sum(e[:6]),sum(e[6:]))
    if bi in pieces:pieces[bi][e]=pieces[bi].get(e,OK.zero())+sign*co*co2
 pieces={bi:OR({e:c for e,c in dd.items() if c}) for bi,dd in pieces.items()}
 cofs=[cofactor(v) for v in hm['items'][n]['cofactors']]
 # The normalized binomials pulled back are the literal L*RR generators.
 bins=[cofactor(v) for v in hm['binomial_generators']]
 ck('fresh complete H32 source module '+str(n),sum(cofs[j]*bins[j] for j in range(3))==pieces[3,2])
 ck('literal L binomials '+str(n),vector(OR,bins)==ell)
 zero=pull(pieces[4,1])+sum(pull(cofs[j])*ee[j] for j in range(3));delta=pull(pieces[1,3]);rows.append(row(zero,delta));data=dict(coefficients=rows[-1],source_faces=pieces,cofactors=cofs,checks=checks[begin_checks:],fingerprints=fingerprints);savedH.append(data);save(data,str(checkpoint),compress=False)
 print('fresh critical fullH',n,round(time.monotonic()-st,2),flush=True)
if partial:
 print('FRESH CRITICAL ROW CHUNK COMPLETE',bounds,flush=True);sys.exit(0)
# Entire even coefficient module witnesses, multiplied against fresh quartics.
ev=load(str(research/'round10/independent_delta_subcritical_even_full_source.sobj'));ek=ev['ring'].base_ring();ep=ek.ring().hom([OK.gen(0),(OK.gen(0)-OK.gen(1))/(OK.gen(0)*(OK.gen(1)-1)),(OK.gen(0)-OK.gen(1))/(OK.gen(0)*(OK.gen(1)-1))*OK.gen(2),(OK.gen(0)-OK.gen(1))/(OK.gen(1)*(OK.gen(0)-1))*OK.gen(3)],OK)
def ec(c):c=ek(c);return ep(c.numerator())/ep(c.denominator())
def emap(v):return OR({tuple(map(int,e))+(0,)*(9-len(e)):ec(c) for e,c in v.dict().items()})
evrows=[]
for item in ev['items']:
 idx=int(item['Q_index']);jj=int(item['square_component']);source=joint['quartic_source_rows'][idx];unit=ec(item['square_unit']);A=[emap(v) for v in item['A_p_quadratic']]
 ck('fresh even22 square '+str(idx),part(source,(2,2))==unit*ell[jj]**2)
 ck('fresh even31 entiremodule '+str(idx),part(source,(3,1))==sum(A[j]*ell[j] for j in range(3)))
 evrows.append(cv(unit)*ee[jj]**2+sum(pull(A[j])*ee[j] for j in range(3))+pull(part(source,(4,0))))
M=matrix(K,rows);save(dict(field=K,e_ring=P,p_RR_boundary=vals,matrix=M,H_items=savedH,even_sources=evrows),str(out/'fresh_critical_linear_sources.sobj'),compress=False)
old=load(str(research/'round10/independent_delta_critical_linear_sources.sobj'));oldK=old['field'];mp=oldK.ring().hom(list(K.gens()),K)
def oc(c):c=oldK(c);return mp(c.numerator())/mp(c.denominator())
for i in range(22):
 for j in range(5):ck('all22x5 critical source coefficients '+str((i,j)),M[i,j]==oc(old['matrix'][i,j]))
for i,v in enumerate(old['even_sources']):ck('all even coefficients '+str(i),evrows[i]==P({tuple(map(int,e)):oc(c) for e,c in v.dict().items()}))
record=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,full22rows=True,all3even=True,scope='Entire critical original source matrix freshly reconstructed from current Section2-derived R/Q/H with every imported module cofactor multiplied exactly; no old source matrix or stored zero used in generation.')
(out/'fresh_critical_linear_sources.json').write_text(json.dumps(record,indent=2));print('PASS',len(checks),record['seconds'],flush=True)
