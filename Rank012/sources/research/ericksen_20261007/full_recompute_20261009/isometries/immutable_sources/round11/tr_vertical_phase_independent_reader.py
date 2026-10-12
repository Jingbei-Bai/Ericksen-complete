"""Independent source and protected-division audit for the pure T-R phase fiber."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent;checks=[]
def ck(label,yes):
 assert yes,label
 checks.append(label)
d=load(str(OUT/'tr_pure_rotation_first_order_fields.sobj'));S=d['P'];K=S.fraction_field()
v,j,A,T,B,L=map(K,S.gens());cols=d['X'];bs=load(str(OUT/'tr_pure_rotation_B_fields.sobj'))['Bz']
for index,col in enumerate(cols):
 ck('Bz original first-order source '+str(index),K(bs[index])==2*B*col[1]+2*B*L/(v+L*L)*col[3])
n=K(bs[0]).numerator();ck('constant source numerator L degree at most three',n.degree(S.gen(5))<=3)
R=PolynomialRing(QQ,names=('v','j','I1','I2','B'));F=R.fraction_field()
vv,jj,i1,i2,bb=map(F,R.gens());pp=vv**3-i1*vv**2+i2*vv-jj;aa=jj/vv-pp/(vv*bb)
theta=vv*aa*(i2-aa-jj*(1+bb)/vv)-jj*vv*vv
actual=[]
for power in [0,3]:
 even=F(0);odd=F(0)
 for ex,co in n.coefficient({S.gen(5):power}).dict().items():
  z=co*vv**ex[0]*jj**ex[1]*aa**ex[2]*theta**(ex[3]//2)*bb**ex[4]
  if ex[3]%2:odd+=z
  else:even+=z
 actual.append((even,odd))
ck('L zero coefficient exactly odd T',actual[0][0]==0)
ck('L cubic coefficient exactly even T',actual[1][1]==0)
z0=F(actual[0][1]/(aa**3*theta));z3=F(actual[1][0]/(aa**3*theta))
ck('z0 exact displayed identity',z0==-4*vv**2*(5*vv**3-3*i1*vv**2+i2*vv+jj*(bb+1)))
for lab,z in [('z0',z0),('z3',z3)]:
 nn=z.numerator();dd=z.denominator();gg=nn.gcd(dd);nn=nn//gg;dd=dd//gg
 ck(lab+' remaining denominator only B monomial',len(dd.dict())==1 and all(tuple(ex[:4])==(0,0,0,0) for ex in dd.dict()))
 if lab=='z0':n0=nn
 else:n3=nn
res=n0.resultant(n3,R.gen(4))
expected=-256*jj**2*vv**9*(-25*vv**3+25*i1*vv**2-(6*i1*i1+7*i2)*vv+3*i1*i2-2*jj)*(5*vv**3-3*i1*vv**2+i2*vv+jj)
ck('entire B resultant exact',F(res)==expected)
ck('highest v coefficient universally nonzero',res.coefficient({R.gen(0):15})==32000*R.gen(1)**2)
report={'status':'PASS','checks':len(checks),'labels':checks,
 'all_divisions':'T to remove the odd part, A^3*T^2, and B powers; all nonzero on actual chart. No pressure or material discriminant divided.',
 'v_exception':'A fixed nonzero degree-six polynomial after protected v,j factors, for every positive spectrum.'}
(OUT/'tr_vertical_phase_independent_reader.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2),flush=True)
