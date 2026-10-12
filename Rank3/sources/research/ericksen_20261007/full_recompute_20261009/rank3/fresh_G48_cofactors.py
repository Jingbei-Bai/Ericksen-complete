"""Regenerate G48 using literal polynomial 3x3 cofactors and Singular resultant.
This is an alternative algorithm to fraction-field matrix inversion.
"""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;start=time.monotonic();checks=[]
src=out.parent.parent/'rank3_simplification_20261009/original_moment_definitions.json'
C=PolynomialRing(QQ,names=('H','Z','S'));H,Z,S=C.gens();K=C.fraction_field();A=PolynomialRing(K,'w');w=A.gen()
env=dict(H=A(H),Z=A(Z),S=A(S),w=w)
data=json.loads(src.read_text());raw=matrix(A,[[sage_eval(e,locals=env) for e in row] for row in data['general_curve_trilinear_canonical_matrix_plain']])
P=PolynomialRing(C,'w');wp=P.gen();units=[];rows=[]
def ck(name,value):
 if not value:raise AssertionError(name)
 checks.append(name)
def pp(p):
 p=A(p);den=lcm([c.denominator() for c in p.list()]);co=gcd([C(c*den) for c in p.list()]);return P([C(c*den)//co for c in p.list()])
for i in range(4):
 den=lcm([c.denominator() for f in raw.row(i) for c in f.list()]);co=gcd([C(c*den) for f in raw.row(i) for c in f.list()]);unit=K(co/den);units.append(unit)
 row=[P([C(c/unit) for c in f.list()]) for f in raw.row(i)];rows.append(row)
 ck('row unit multiplication '+str(i),all(A(row[j])*unit==raw[i,j] for j in range(4)))
N=matrix(P,rows)
def det3(B):return B[0,0]*(B[1,1]*B[2,2]-B[1,2]*B[2,1])-B[0,1]*(B[1,0]*B[2,2]-B[1,2]*B[2,0])+B[0,2]*(B[1,0]*B[2,1]-B[1,1]*B[2,0])
adj=matrix(P,4,4)
for i in range(4):
 for j in range(4):adj[i,j]=(-1)**(i+j)*det3(N.matrix_from_rows_and_columns([x for x in range(4) if x!=j],[x for x in range(4) if x!=i]))
 print('cofactor row',i,'seconds',time.monotonic()-start,flush=True)
det=sum(N[0,j]*adj[j,0] for j in range(4));ck('cofactor adjugate multiplication',N*adj==det*identity_matrix(P,4));D=pp(det)
L=H*S-Z;Q=(2*H*Z*wp+L)*(2*Z*(1+H-Z)*wp-(S-1)*(L+S));RR=PolynomialRing(QQ,names=('H','Z','S','w'));HH,ZZ,SS,WW=RR.gens()
def globalp(p):return sum(RR({(int(e[0]),int(e[1]),int(e[2]),i):c for e,c in cc.dict().items()}) for i,cc in enumerate(p.list()))
DG=globalp(D);res=[];quartics=[]
for j in [0,1]:
 col=adj.column(j);f=pp(Q*col[1]*col[2]-col[0]*col[3]);quartics.append(f);FG=globalp(f)
 ck('fixed quartic degree '+str(j),f.degree()==4)
 print('resultant start',j,'degree',D.degree(),f.degree(),'seconds',time.monotonic()-start,flush=True)
 rr=DG.resultant(FG,WW);rr=C({tuple(e[:3]):c for e,c in rr.dict().items()});res.append(rr)
 print('resultant end',j,'terms',len(rr.dict()),'seconds',time.monotonic()-start,flush=True)
B=H**3*Z*S*S-H*Z**3-H*Z*S*S-2*H*H*S+2*Z*Z*S+H*Z;guard=C(H*Z*(H-Z)*(Z-1)*(H+1)*(1+H-Z)*(S-1)*(H*S-Z+1)*((H+1)*S-Z)*S*(H*S-Z))
gc=gcd(*res);ck('B divides gcd',gc.quo_rem(B)[1]==0);G=gc//B
for factor,power in guard.factor():
 while G.quo_rem(factor)[1]==0:G=G//factor
G=C(G/G.leading_coefficient());ck('G48 total degree',G.total_degree()==48);ck('G S degree14',G.degree(S)==14)
Fd=load(str(out/'fresh_F38.sobj'));F=C(Fd['F']);contact=F.derivative(S)*(G.derivative(H)+S*G.derivative(Z))-G.derivative(S)*(F.derivative(H)+S*F.derivative(Z));ck('contact degree25',contact.degree(S)==25)
save(dict(ring=C,F=F,G=G,Ccontact=contact,sources=[F,G,contact],D=D,resultants=res,quartics=quartics,guard=guard,row_units=units,normalized_adj=adj),str(out/'fresh_FG_sources.sobj'),compress=False)
(out/'fresh_FG_sparse.json').write_text(json.dumps(dict(variables=['H','Z','S'],F=[[list(map(int,e)),str(c)] for e,c in F.dict().items()],G=[[list(map(int,e)),str(c)] for e,c in G.dict().items()],Ccontact=[[list(map(int,e)),str(c)] for e,c in contact.dict().items()]),indent=2))
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-start,source='Explicit manuscript matrices; row units retained; direct polynomial cofactors; fixed-degree2,4 Singular resultant',Fdegree=int(F.total_degree()),Gdegree=int(G.total_degree()),Fterms=len(F.dict()),Gterms=len(G.dict()),no_saved_resultant_or_factor_input=True)
(out/'fresh_G48_cofactors.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
