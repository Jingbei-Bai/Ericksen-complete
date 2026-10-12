"""Rebuild all 21 rank-three quartic sources from fresh Section 2 fields.

The only mathematical inputs are the freshly rebuilt 17-state system and the
three explicit RR coordinate formulas printed in Section 5. Cubic elimination,
the combined u reconstruction, all quartic rows, and all three actual state
fields are recomputed.
"""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;start=time.monotonic();checks=[]
fresh=load(str(out.parent/'foundations/state17_fresh.sobj'))
P17=fresh['ring'];Ko=P17.base_ring()
CK=PolynomialRing(QQ,names=('a','b','h','z'));K=CK.fraction_field();a,b,h,z=K.gens()
R=PolynomialRing(K,names=('p12','p13','p21','p23','p31','p32','rr1','rr2','rr3'))
px=R.gens();p12,p13,p21,p23,p31,p32,R1,R2,R3=px
T=PolynomialRing(R,names=('t0','t4','t5','u0','u1'));t0,t4,t5,u0,u1=T.gens()
def ck(name,condition):
 if not condition:raise AssertionError(name)
 checks.append(name)
def prog(name,**kw):print(name,round(time.monotonic()-start,2),kw,flush=True)
oldco=Ko.ring().hom([a,b],K)
def cv(c):c=Ko(c);return oldco(c.numerator())/oldco(c.denominator())
def evaluate(f,images,target):
 outp=target.zero()
 for e,c in f.dict().items():
  term=target(cv(c))
  for nn,img in zip(e,images):
   if nn:term*=img**int(nn)
  outp+=term
 return outp
r=(a-b)/(a*(b-1));q=1+r;H=r*h;Z=r*z/q
t1=-p12*(p13+p31)+(1-H*Z/(H-Z))*p23-H*Z/(H-Z)*p32-H*Z/(2*(H-Z))*R1
t2=-2*H/(r*(1+H))*p12+2*H*H/(r*(1+H))*p21-H/(r*(1+H))*R3
t3=p21*(p23+p32)+(Z-H+H*Z)/(r*(1-Z))*p13+Z/(r*(1-Z))*p31+Z/(2*r*(1-Z))*R2
images=list(map(T,[1,p12,p13,p21,h,p23,p31,p32,z,t0,t1,t2,t3,t4,t5,u0,u1]))
cs=[evaluate(f,images,T) for f in fresh['C3']]
mat=matrix(K,[[K(f.derivative(v).constant_coefficient()) for v in (t0,t4,t5)] for f in cs])
ck('cubic coefficient matrix contains no state variables',all(cs[i].derivative(v)==mat[i,j] for i in range(3) for j,v in enumerate((t0,t4,t5))))
rhs=vector(R,[-f(t0=0,t4=0,t5=0,u0=0,u1=0) for f in cs]);ts=mat.inverse()*rhs
for index,value in zip([9,13,14],ts):images[index]=T(value)
ck('all cubic sources vanish on reconstructed graph',all(evaluate(f,images,T)==0 for f in fresh['C3']))
prog('fresh cubic reconstruction')
qq=[evaluate(f,images,T) for f in fresh['Q4']]
M=matrix(R,[[R(f.derivative(v).constant_coefficient()) for v in (u0,u1)] for f in qq])
A=a*(b-1);C=b*(a-1);D=a-b
B0=M.matrix_from_rows([0,3]);B1=M.matrix_from_rows([12,16])
comb=A/2*B0.adjugate()*B0+2*A*A/D*B1.adjugate()*B1
ck('joint combined u left inverse',comb==7*C*h*identity_matrix(R,2))
rhsu=vector(R,[f(u0=0,u1=0,t0=0,t4=0,t5=0) for f in qq])
us=-(A/2*B0.adjugate()*vector(R,[rhsu[i] for i in [0,3]])+2*A*A/D*B1.adjugate()*vector(R,[rhsu[i] for i in [12,16]]))/(7*C*h)
images[15]=T(us[0]);images[16]=T(us[1]);images=[R(v) for v in images]
rows=[evaluate(f,images,R) for f in fresh['Q4']]
ck('all quartics at most quadratic RR',all(all(sum(e[6:])<=2 for e in f.dict()) for f in rows))
ck('RR Jacobian protected nonzero',matrix(K,[[K(v.derivative(x)) for x in px[6:]] for v in [t1,t2,t3]]).det()!=0)
prog('all 21 quartics independently rebuilt')
save(dict(ring=R,parameters=(a,b,h,z),variables=px,quartic_source_rows=rows,source='Fresh Section 2 state and printed RR reconstruction'),str(out/'fresh_joint_quartics.sobj'),compress=False)
raw=[[evaluate(f,images,R) for f in row] for row in fresh['F']]
weights=[1]*9+[2]*6+[3]*2
norm=[[raw[i][j]-weights[j]*images[j]*raw[i][0] for j in range(17)] for i in range(3)]
idx=[1,2,3,5,6,7];tf=[t1,t2,t3];df=matrix(K,[[K(t.derivative(v)) for v in px[6:]] for t in tf]);di=df.inverse()
fields=[];bases=[]
for i in range(3):
 off=vector(R,[norm[i][j] for j in idx]);base=vector(R,[norm[i][4],norm[i][8]])
 chain=vector(R,[sum(t.derivative(px[j])*off[j] for j in range(6))+R({e:c.derivative(h) for e,c in t.dict().items()})*base[0]+R({e:c.derivative(z) for e,c in t.dict().items()})*base[1] for t in tf])
 vf=di*(vector(R,[norm[i][j] for j in [10,11,12]])-chain)
 ck('fresh joint field chain '+str(i),df*vf+chain==vector(R,[norm[i][j] for j in [10,11,12]]))
 fields.append(list(off)+list(vf));bases.append(list(base));prog('fresh state direction',direction=i)
save(dict(ring=R,parameters=(a,b,h,z),variables=px,fields=matrix(R,fields),base_fields=matrix(R,bases),quartic_source_rows=rows,source='Fresh Section 2 state and printed RR reconstruction'),str(out/'fresh_joint_fields.sobj'),compress=False)
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-start,source_generation='Fresh state17, printed RR formulas; independent cubic and u solves',historical_inputs='None',no_material_specialization=True)
(out/'fresh_joint_sources.json').write_text(json.dumps(report,indent=2));prog('PASS',count=len(checks))
