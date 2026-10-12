"""All necessary full sources on the entire S=N=0 curve family.

Z is an arbitrary nonzero constant with all inherited guards retained;
H remains a variable base coordinate. No spectrum is specialized.
"""
from pathlib import Path
import time,json,sys
from sage.all import *
out=Path(__file__).resolve().parent;start=time.monotonic()
full=load(str(out/'general_curve_dimensionless_full.sobj'))
extra=load(str(out/'general_curve_extra_quartic_cubics.sobj'))
q4=load(str(out/'general_curve_remaining_quartics.sobj'))
proportional='--proportional' in sys.argv
stem='general_curve_proportional_Z_full_sources' if proportional else 'general_curve_constant_Z_full_sources'
oldR=full['ring'];C=PolynomialRing(QQ,names=('a','r','H','c' if proportional else 'Z'));K=C.fraction_field();a,r,H,curve_parameter=K.gens()
Z=curve_parameter*H if proportional else curve_parameter
S=curve_parameter if proportional else K.zero()
R=PolynomialRing(K,names=('x1','y1','x2','y2','x3','y3'));x1,y1,x2,y2,x3,y3=R.gens()
# Original, unscaled opposite offdiagonal entries; no square root used.
im=[x1,-2*(H*y1+Z*x1)/(1-S),x2,-2*H*(y2+Z*x2)/(1+H*S-Z),x3,2*Z*(y3-H*x3)/((1+H)*S-Z)]
T=PolynomialRing(K,names=oldR.variable_names());ph=T.hom(im,R)
def cc(c):return K(c.numerator()(a,r,H,Z,S,0))/K(c.denominator()(a,r,H,Z,S,0))
def ev(f):return ph(T({m:cc(c) for m,c in f.dict().items()}))
sources=full['rows']+extra['rows']+q4['rows'];rows=[];primitive=[];scalings=[];records=[]
names=['mixed01','mixed02','mixed12','Q0','S00','S11','S22']+['Q'+str(t[0])+'_cubic' for t in extra['row_specs']]+['Q1','Q5','Q6']
for name,p in zip(names,sources):
 f=ev(p);rows.append(f);den=lcm([c.denominator()/c.denominator().leading_coefficient() for c in f.coefficients()]);den/=den.leading_coefficient();num=[C(den*c) for c in f.coefficients()];con=gcd(num)
 prim=R({m:C(c*den)//con for m,c in f.dict().items()});primitive.append(prim);scalings.append((con,den))
 hist={}
 for mon in f.dict():hist[sum(mon)]=hist.get(sum(mon),0)+1
 rec={'source':name,'terms':int(f.number_of_terms()),'degrees':hist,'content':str(con.factor()),'denominator':str(den.factor())}
 records.append(rec);print(json.dumps(rec),flush=True)
payload={'ring':R,'rows':rows,'primitive_rows':primitive,'scalings':scalings,'names':names,'images':im,'records':records,'parameter_images':(a,r,H,Z,S,0),'proportional':proportional,'scope':'Full necessary sources for the entire constant-Z or Z=cH noncharacteristic actual ratio curve family. Curve parameter and full a,r spectrum retained; no dimension conclusion.'}
save(payload,str(out/(stem+'.sobj')))
(out/(stem+'.json')).write_text(json.dumps({'records':records,'elapsed':time.monotonic()-start},indent=2)+'\n')
