"""Fresh source matrices for constant-column cover and two Cramer certificates.
Certificate multipliers are witnesses; every multiplication and their source
provenance are recomputed against the freshly derived manuscript matrix.
"""
from sage.all import *
from itertools import combinations
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;research=out.parent.parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
def progress(n):print(n,round(time.monotonic()-st,2),flush=True)
d=load(str(out/'fresh_source_to_moments.sobj'));Ms=d['printed_matrix'];A=Ms.base_ring();K=A.base_ring();C=K.ring();H,Z,S=K.gens();w=A.gen();L=H*S-Z;cap=1+H-Z
M=matrix(A,4,4,[Ms[i,j](w+L/(2*H*Z)) for i in range(4) for j in range(4)])
Q=(2*H*Z*w+L)*(2*Z*cap*w-(S-1)*(L+S));aq=Q[2];bq=Q[1]
P=matrix(A,[[1,H*(w+bq/aq),-w,H/aq],[0,1,0,0],[0,-H,1,0],[0,-aq*w,0,1]]);NN=M*P
for i in range(4):
 for j in range(4):ck('fresh column degree '+str((i,j)),NN[i,j].degree()<=[1,0,1,0][j])
atoms=[C(v) for v in [H,Z,H-Z,Z-1,H+1,cap,S-1,L+1,L+S,S,L]];guard=prod(atoms)
def strip(poly):
 p=C(poly)
 for atom in atoms:
  while p and p.quo_rem(atom)[1]==0:p=p//atom
 return C(p/p.leading_coefficient()) if p else p
def only_guard(c):
 c=K(c)
 return bool(c) and all(guard.quo_rem(g)[1]==0 for part in [c.numerator(),c.denominator()] for g,e in part.factor())
certfile=research/'round10/general_curve_mixed_short_unit_integer.json';cert=json.loads(certfile.read_text())
def pol(data):return C({tuple(map(int,e)):QQ(c) for e,c in data})
sources=[pol(f) for f in cert['sources']];mults=[pol(f) for f in cert['multipliers']];DH=pol(cert['H_multiplier']);tar=pol(cert['target'])
ck('portable exact integer identity',sum(a*p for a,p in zip(mults,sources))==DH*tar)
ck('H denominator nonzero H only',DH!=0 and DH.degree(C.gen(1))==DH.degree(C.gen(2))==0)
ck('exact target S2Sminus12',tar==C(S*S*(S-1)**2))
fd=load(str(out/'fresh_F38.sobj'));F=C(fd['F']);deltas=[C(p) for p in fd['deltas']]
for i,(delta,sr) in enumerate(zip(deltas,sources)):
 ratio=K(delta)/(F*sr);ck('fresh determinant certificate source '+str(i),only_guard(ratio))
constant=matrix(K,[[NN[i,j][0] for j in [1,3]] for i in range(4)])
for i,(rows,residual) in enumerate(zip(combinations(range(4),2),[0,1,3,2,4,5])):
 minor=constant.matrix_from_rows(rows).det();ck('fresh constant-column source '+str(i),only_guard(minor/sources[residual]))
progress('source-to-integer-certificate and six pivot cover')
linear=[]
for rows in combinations(range(4),3):
 for col in [0,2]:
  p=NN.matrix_from_rows_and_columns(rows,[1,3,col]).det();ck('lowrank pencil linear',p.degree()<=1)
  den=lcm([v.denominator()/v.denominator().leading_coefficient() for v in p.list()]);con=gcd([C(v*den) for v in p.list()]);prot=C.one()
  for atom in atoms:
   while con and con.quo_rem(atom)[1]==0:con=con//atom;prot*=atom
  pp=A([K(C(v*den)//prot) for v in p.list()]);scale=pp[0].numerator().leading_coefficient() if pp[0] else pp[1].numerator().leading_coefficient();pp/=scale;prot*=scale
  ck('full lowrank source identity',den*p==prot*pp and only_guard(K(prot/den)));linear.append(pp)
BB=matrix(C,[[C(f[i]) for i in range(2)] for f in linear]);parameters=[strip(BB.matrix_from_rows(rows).det()) for rows in combinations(range(8),2)]
arch=load(str(research/'round10/general_curve_canonical_ranktwo_linear.sobj'))
for i,(f,g) in enumerate(zip(parameters,arch['parameter_sources'])):ck('fresh ranktwo source equals Cramer input '+str(i),f==C(g))
zero=load(str(research/'round10/general_curve_zero_det_coefficients_input.sobj'));newD=load(str(out/'fresh_FG_sources.sobj'))['D']
for i,item in enumerate(zero['items']):
 new=C(newD[i]);expected=C(item['coefficient']);ck('fresh determinant coefficient matches input '+str(i),K(new/expected) in QQ and new!=0)
 ck('fresh determinant stripped source '+str(i),only_guard(K(new)/C(item['source'])))
save(dict(constant_column_matrix=constant,lowrank_linear_sources=linear,lowrank_parameter_sources=parameters,unit_sources=sources,unit_multipliers=mults,H_multiplier=DH,target=tar),str(out/'fresh_rank_matrix_certificates.sobj'),compress=False)
rep=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st,lowrank_sources=28,zero_det_sources=3,constant_pivot_charts=6,unit_witness_replayed=True,
 scope='Fresh actual-source matrices equal all original coefficient sources of the ranktwo/zero-det Cramer certificates. Complete integer-node/degree replay is in separate isolated exact-source readers; not a new rational solution computation.')
(out/'fresh_rank_matrix_certificates.json').write_text(json.dumps(rep,indent=2));progress('PASS')
