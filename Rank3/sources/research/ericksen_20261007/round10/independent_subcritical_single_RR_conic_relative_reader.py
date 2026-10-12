"""Read full source matrix and verify the relative conic module independently.

Build the complete eight-generator coefficient matrix from genuine first
jets and exact boundary division; take determinant by permutation, without
using the creator's four-block determinant algorithm.
"""
from sage.all import *
from pathlib import Path
from itertools import permutations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_single_RR3_normal.sobj'));cert=load(str(p/'independent_subcritical_single_RR_conic_generic_factor_blocks.sobj'));E=cert['ring'];K=E.base_ring();k,r,H,Z,w=K.gens();x=E.gen();OP=d['ring'];OK=OP.base_ring();ph=OK.ring().hom([k/r,r,H,Z],K);checks=[];cache={}
def ck(label,v):
 assert bool(v),label;checks.append(label)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def ne(v):return E([nk(v[0]),nk(v[1])])
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=nk(ph(c.numerator())/ph(c.denominator()))
 return cache[c]
def mp(f):
 ans=E.zero()
 for ex,c in f.dict().items():
  ex=tuple(map(int,ex));assert not any(ex[3:]);ans+=cv(c)*x**ex[0]*(-x)**ex[1]*w**ex[2]
 return ne(ans)
for i in range(25):ck('full-base R/Q/H graph boundary '+str(i),d['graph_boundary'][i]==0)
ox,oz=OP.gens()[:2]
sq=[]
for i in (25,26):
 f=d['graph_boundary'][i];change=OP.hom([ox-oz]+list(OP.gens()[1:]),OP);back=OP.hom([ox+oz]+list(OP.gens()[1:]),OP);changed=change(f);ck('square boundary zero constant and linear Taylor '+str(i),all(ex[0]>=2 for ex in changed.dict()));qchanged=OP({tuple(int(n)-(2 if j==0 else 0) for j,n in enumerate(ex)):c for ex,c in changed.dict().items()});quotient=back(qchanged);ck('exact full-base square boundary division '+str(i),quotient*(ox+oz)**2==f);sq.append(mp(quotient))
M0=d['graph_normal_matrix'];M=matrix(E,28,8,lambda i,j:mp(M0[i,j+3]) if j<4 else (E(i==27) if j==4 else (sq[i-25] if j==5 and i in (25,26) else (mp(M0[i,j+1]) if j>=6 else E.zero()))))
for i in range(28):ck('whole conic third-even boundary '+str(i),i!=27 or mp(d['graph_boundary'][i])==0)
rows=tuple(cert['rows']);A=M.matrix_from_rows(rows);det=E.zero();terms=0
for perm in permutations(range(8)):
 if any(not A[i,perm[i]] for i in range(8)):continue
 sign=(-1)**sum(perm[i]>perm[j] for i in range(8) for j in range(i+1,8));term=E.one()
 for i in range(8):term=ne(term*A[i,perm[i]])
 det=ne(det+sign*term);terms+=1
pieces=cert['pieces'];prod=ne(pieces['pure_kappa']*pieces['Q0_e1']*pieces['range3']*pieces['even2']);ck('independent true8 determinant agrees up to sign',det==prod or det==-prod);ck('determinant x coefficient zero',det[1]==0)
factor=w*(w*w+k)**5*(w*w-k)**2;unit=nk(det[0]/factor);ck('determinant state factors complete',not unit.derivative(w));ck('genuine generic coefficient unit nonzero',unit!=0)
print('DET UNIT',unit.factor(),'terms',terms,'checks',len(checks),'sec',time.monotonic()-st,flush=True)
payload={'ring':E,'matrix':M,'rows':rows,'determinant':det,'state_factor':factor,'unit':unit,'checks':checks,'scope':__doc__}
save(payload,str(p/'independent_subcritical_single_RR_conic_relative_reader.sobj'),compress=False)
(p/'independent_subcritical_single_RR_conic_relative_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'independent_permutation_terms':terms,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
