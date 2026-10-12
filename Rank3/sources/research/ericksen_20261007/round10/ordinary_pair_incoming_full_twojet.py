"""Literal three complete pair sources on the ordinary incoming family.

The entire polynomial tau^4*S(p/tau,RR/tau) is reconstructed before
extracting its boundary and first/second normal coefficients. No source
vanishing, material specialization, or range correction is presumed.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations_with_replacement
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'independent_ratio_p_boundary_global_pair_syzygies.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+600*1024**2,)*2)
old=d['coefficient_field'];K=PolynomialRing(QQ,names=('op_a','op_r','op_H','op_Z')).fraction_field();a,r,H,Z=K.gens();q=1+r;m=1+a*r
ph=old.ring().hom([a,a*q/m,H/r,q*Z/r],K);cache={}
P=PolynomialRing(K,names=('op_l','op_u','op_v','op_n0','op_n1','op_n2','op_n3','op_eta','op_tau'))
l,u,v,n0,n1,n2,n3,eta,tau=P.gens();lambda0=-2*H*(Z+1)/(Z*(H+1));lambda1=(Z-H)/(Z*(H+1));lam=lambda0+lambda1*v
U=PolynomialRing(K,names=('op_l','op_u','op_v'));uv=U.gens()
ne=[(0,)*6]
for j in range(6):
 e=[0]*6;e[j]=1;ne.append(tuple(e))
for i,j in combinations_with_replacement(range(6),2):
 e=[0]*6;e[i]+=1;e[j]+=1;ne.append(tuple(e))
def ck(label,x):
 assert bool(x),label;checks.append(label)
def cv(c):
 c=old(c)
 if c not in cache:
  num=ph(c.numerator());den=ph(c.denominator());assert den
  cache[c]=num/den
 return cache[c]
items=[];mat=[]
for it in d['items']:
 full=it['full'];terms={};src={}
 for ex,oc in full.dict().items():
  e=tuple(map(int,ex));dd=sum(e);ck('nonnegative ordinary homogenization '+str((it['pair'],e)),dd<=4)
  c=cv(oc);src[e]=c
  for je in range(e[6]+1):
   for jv in range(e[6]-je+1):
    coeff=c*binomial(e[6],je)*binomial(e[6]-je,jv)*lambda0**(e[6]-je-jv)*lambda1**jv
    ee=(e[5]+e[6]-je,e[7],e[8]+jv)+e[1:5]+(je,4-dd)
    terms[ee]=terms.get(ee,K.zero())+coeff
 f=P(terms);coeffs={nn:{} for nn in ne}
 for ee,c in f.dict().items():
  e=tuple(map(int,ee));nn=e[3:]
  if nn in coeffs:coeffs[nn][e[:3]]=c
 vv=[U(coeffs[nn]) for nn in ne];mat.append(vv)
 # Independent direct polynomial substitution, without the binomial loop.
 imgs=(P.one(),n0,n1,n2,n3,l,lam*l+eta,u,v)
 direct=sum((c*tau**(4-sum(e))*prod(imgs[j]**e[j] for j in range(9)) for e,c in src.items()),P.zero())
 ck('full ordinary source reconstruction '+str(it['pair']),f==direct)
 boundary=vv[0]
 print('PAIR',it['pair'],'FULL',len(f.dict()),'BOUNDARY',str(boundary.factor()) if boundary else '0','FIRST',[len(z.dict()) for z in vv[1:7]],'SECONDS',[len(z.dict()) for z in vv[7:]],'sec',time.monotonic()-st,flush=True)
 items.append({'pair':it['pair'],'full_homogeneous_pullback':f,'boundary':boundary,'source_full_mapped':src,'source_combination':it['full_source_combination'],'jet_coefficients':vv})
 save({'coefficient_field':K,'ring':P,'jet_ring':U,'lambda':lam,'normal_positions':('p13','p21','p23','p31','eta','tau'),'normal_exponents':ne,'items':items,'coefficient_matrix':matrix(U,mat),'checks':checks,'scope':__doc__},str(p/'ordinary_pair_incoming_full_twojet.sobj'),compress=False)
rec={'status':'FULL ORIGINAL SOURCES AND JETS REPLAYED','checks':len(checks),'seconds':time.monotonic()-st,'boundary_zero':[not bool(it['boundary']) for it in items]}
(p/'ordinary_pair_incoming_full_twojet.json').write_text(json.dumps(rec,indent=2)+'\n');print('DONE',rec,flush=True)
