"""Simple graph on the persistent B wall: no B_S inverse is required."""
from pathlib import Path
import json,time
from sage.all import load,save,PolynomialRing,QQ
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'general_curve_B_wall_contact.sobj'))
C=d['ring'];H,Z,S=C.gens();L=H*S-Z;B=H*Z*(H**2-1)*S**2+2*(Z**2-H**2)*S-H*Z*(Z**2-1)
F=H*Z*(B.derivative(H)+S*B.derivative(Z))+S*L*B.derivative(S)
checks={'B_source':B==d['B'],'polynomial_contact_no_inverse':F==d['simple_contact_numerator']}
P=PolynomialRing(QQ,names=('H','Z'));h,z=P.gens();R=PolynomialRing(P.fraction_field(),'ss');ss=R.gen()
def cv(f):return sum((c*h**int(m[0])*z**int(m[1])*ss**int(m[2]) for m,c in f.dict().items()),R.zero())
res=cv(B).resultant(cv(F));expected=-12*(z-h)*(z-1)*(z+1)*(h-1)*(h+1)*(h+z)*z**5*h**5*(h**2+z**2+1)
checks['complete_resultant_recomputed']=res==expected==d['simple_contact_resultant']
def contact(f):return f.derivative(H)+S*f.derivative(Z)
checks['Zplusone_contact']=contact(Z+1)==S
checks['ZplusH_contact']=contact(Z+H)==S+1
checks['ZequalsminusH_contact_forces_Lzero']=L.subs({Z:-H,S:-1})==0
disc=(Z**2-H**2)**2+H**2*Z**2*(H**2-1)*(Z**2-1)
sos=((H**2*Z**2-H**2)**2+(H**2-Z**2)**2+(Z**2-H**2*Z**2)**2)/2
checks['discriminant_source_and_SOS']=disc==sos==d['discriminant_quarter']
checks['quadratic_discriminant_identity_no_A_inverse']=B.derivative(S)**2-4*H*Z*(H**2-1)*B==4*disc
checks['B_and_BS_at_unique_allowed_real_zero']=B.subs({H:1,Z:-1})==0 and B.derivative(S).subs({H:1,Z:-1})==0
assert all(checks.values())
scope='Simple graph cannot persist on B=0 along a real H-submersion interval. B_S is nonzero on an open dense subset of any real actual B=0 interval, by the exact discriminant SOS. Other directions on B=0 are not excluded.'
save({'checks':checks,'B':B,'F':F,'resultant':res,'scope':scope},str(out/'independent_B_wall_simple_contact.sobj'))
(out/'independent_B_wall_simple_contact.json').write_text(json.dumps({'checks_passed':len(checks),'scope':scope,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
print({'checks_passed':len(checks),'seconds':time.monotonic()-start},flush=True)
