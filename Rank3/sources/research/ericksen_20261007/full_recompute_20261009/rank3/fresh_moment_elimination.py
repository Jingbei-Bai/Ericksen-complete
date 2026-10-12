"""Regenerate moment minors, F38 and G48 from explicit printed matrix entries.

The input is a textual definition of the manuscript matrices, not computed
minors, factors, roots or saved PASS records. No material is specialized.
"""
from sage.all import *
from pathlib import Path
from itertools import combinations
import json,time,hashlib,sys
out=Path(__file__).absolute().parent
source=out.parent.parent/'rank3_simplification_20261009'/'original_moment_definitions.json'
start=time.monotonic();checks=[]
def ck(name,value):
 if not value:raise AssertionError(name)
 checks.append(name)
def report(phase,**kwargs):
 data=dict(phase=phase,seconds=time.monotonic()-start,checks=checks,**kwargs)
 (out/'fresh_moment_progress.json').write_text(json.dumps(data,indent=2))
 print(phase,round(data['seconds'],2),kwargs,flush=True)
C=PolynomialRing(QQ,names=('H','Z','S'),order='degrevlex'); H,Z,S=C.gens();K=C.fraction_field()
A=PolynomialRing(K,'w');w=A.gen();env={'H':A(H),'Z':A(Z),'S':A(S),'w':w}
data=json.loads(source.read_text())
def parse_matrix(key):return matrix(A,[[sage_eval(e,locals=env) for e in row] for row in data[key]])
Ms=parse_matrix('general_curve_mixed_Q_matrix_plain');M=parse_matrix('general_curve_trilinear_canonical_matrix_plain')
guard=C(H*Z*(H-Z)*(Z-1)*(H+1)*(1+H-Z)*(S-1)*(H*S-Z+1)*((H+1)*S-Z)*S*(H*S-Z))
guards=[f for f,e in guard.factor()]
def strip(f):
 f=C(f)
 for g in guards:
  while f and f.quo_rem(g)[1]==0:f=f//g
 return f/f.leading_coefficient() if f else f
def primitive(p):
 den=lcm([v.denominator() for v in p.list()]);nums=[C(v*den) for v in p.list()];con=gcd(nums)
 ans=A([K(v//con) for v in nums])
 ck('primitive multiplication '+str(len(checks)),con*ans==den*p)
 return ans,con,den
L=H*S-Z;shift=K(L/(2*H*Z))
for i in range(4):
 for j in range(4):ck('matrix translation '+str((i,j)),Ms[i,j](w+shift)==M[i,j])
report('matrices_loaded',input_sha256=hashlib.sha256(source.read_bytes()).hexdigest())
prims=[];minor_data=[]
for rows in combinations(range(5),4):
 d=Ms.matrix_from_rows(rows).det();pp,con,den=primitive(d)
 ck('maximal minor degree '+str(rows),pp.degree()<=2)
 prims.append(pp);minor_data.append(dict(rows=rows,determinant=d,primitive=pp,content=con,denominator=den))
 report('maximal_minor',rows=list(rows),degree=int(pp.degree()))
MC=matrix(C,[[C(p[j]) for j in range(3)] for p in prims]);deltas=[];gc=None
for rows in combinations(range(5),3):
 d=MC.matrix_from_rows(rows).det();deltas.append(d);gc=d if gc is None else gcd(gc,d)
 report('coefficient_minor',rows=list(rows),gcd_degree=int(gc.total_degree()))
F=C(strip(gc));ck('F total degree',F.total_degree()==38);ck('F S degree',F.degree(S)==11)
save(dict(ring=C,matrix=MC,deltas=deltas,F=F,source_minors=minor_data),str(out/'fresh_F38.sobj'),compress=False)
report('F38_regenerated',terms=len(F.dict()))
if '--full-slow' not in sys.argv:
 report('PASS_F38',status='PASS',Fdegree=int(F.total_degree()),FSdegree=int(F.degree(S)),checks_count=len(checks),next_program='fresh_G48_cofactors.py',scope='All five maximal minors and ten coefficient minors regenerated from textual matrices; upstream source derivation remains separate.')
 (out/'fresh_moment_elimination.json').write_text((out/'fresh_moment_progress.json').read_text())
 sys.exit(0)
D,dcon,dden=primitive(M.det());Adj=M.adjugate();Q=(2*H*Z*w+L)*(2*Z*(1+H-Z)*w-(S-1)*(L+S))
rs=[];fdata=[]
for j in [0,1]:
 col=Adj.column(j);raw=Q*col[1]*col[2]-col[0]*col[3];fp,con,den=primitive(raw)
 ck('quartic degree '+str(j),fp.degree()<=4)
 # Explicit Sylvester determinant ensures fixed degrees (2,4).
 coeffd=[C(D[i]) for i in range(3)];coefff=[C(fp[i]) for i in range(5)]
 Syl=matrix(C,6,6,0)
 for i in range(4):
  for h0 in range(3):Syl[i,i+h0]=coeffd[2-h0]
 for i in range(2):
  for h0 in range(5):Syl[4+i,i+h0]=coefff[4-h0]
 res=Syl.det();rs.append(res);fdata.append(dict(raw=raw,primitive=fp,content=con,denominator=den,resultant=res))
 report('canonical_resultant',column=j,degree=int(res.total_degree()),terms=len(res.dict()))
B=H**3*Z*S*S-H*Z**3-H*Z*S*S-2*H*H*S+2*Z*Z*S+H*Z
gc=gcd(*rs);ck('B divides resultant gcd',gc.quo_rem(B)[1]==0)
G=C(strip(gc//B));ck('G total degree',G.total_degree()==48);ck('G S degree',G.degree(S)==14)
contact=C(F.derivative(S)*(G.derivative(H)+S*G.derivative(Z))-G.derivative(S)*(F.derivative(H)+S*F.derivative(Z)))
ck('contact S degree',contact.degree(S)==25)
save(dict(ring=C,F=F,G=G,Ccontact=contact,sources=[F,G,contact],D=D,adj=Adj,resultants=rs,quartics=fdata,guard=guard),str(out/'fresh_FG_sources.sobj'),compress=False)
(out/'fresh_FG_sparse.json').write_text(json.dumps(dict(variables=['H','Z','S'],F=[[list(map(int,e)),str(c)] for e,c in F.dict().items()],G=[[list(map(int,e)),str(c)] for e,c in G.dict().items()],Ccontact=[[list(map(int,e)),str(c)] for e,c in contact.dict().items()]),indent=2))
report('PASS',status='PASS',Fdegree=int(F.total_degree()),Gdegree=int(G.total_degree()),checks_count=len(checks),no_material_specialization=True,scope='Fresh determinant/minor/adjugate/resultant/gcd calculations from textual manuscript matrices; derivation of matrix rows from physical source fields is separate.')
(out/'fresh_moment_elimination.json').write_text((out/'fresh_moment_progress.json').read_text())
