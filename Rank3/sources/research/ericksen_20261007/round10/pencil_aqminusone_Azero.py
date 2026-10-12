"""Close the leading-coefficient degeneration A=q^2+6q+1=0 on J=0."""
from pathlib import Path
import json
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;d=load(str(out/'pencil_aqminusone_upper_walls.sobj'))
R=PolynomialRing(QQ,'q');q=R.gen();A=q*q+6*q+1
c8=d['coefficients'][8];ans=R.zero()
for mon,val in c8.dict().items():ans+=val*q**int(mon[0])*(QQ(4)/3)**int(mon[1])
rem=ans%A;bez=A.xgcd(rem)
assert bez[0]==1
print('remainder',rem,flush=True);print('gcd',bez[0],flush=True)
save({'wall':A,'coefficient_Hb109':ans,'remainder':rem,'bezout':bez,'T':QQ(4)/3},str(out/'pencil_aqminusone_Azero.sobj'))
(out/'pencil_aqminusone_Azero.json').write_text(json.dumps({'status':'PASS_AZERO_J_WALL','T':'4/3','remainder':str(rem),'gcd':'1','identity_check':bool(bez[1]*A+bez[2]*rem==1)},indent=2)+'\n')
