from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,v):
 if not bool(v):raise AssertionError(n)
 checks.append(n)
d=load(str(p/'independent_subcritical_two_RR_delta_even_elimination.sobj'));P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();t,w,s,e=P.gens();F=d['F38'];even=d['even_sources'];M=load(str(p/'independent_subcritical_two_RR_full_normal.sobj'))['matrix'];mf=M.matrix_from_rows_and_columns([0,4,5,6],[0,1,3,4]).det();ratio=mf/(w*t*F);ck('F is actual original four-minor factor',ratio in K and ratio)
for it in d['necessary_eliminants']:
 j=it['even_row'];v=F.derivative(e)*even[j].subs({e:0})-even[j].derivative(e)*F.subs({e:0});ck('true e resultant without dividing e coefficient',v==it['full_resultant']);den=it['clearing_denominator'];rec=P.zero()
 for hz,f in it['base_coefficients'].items():rec+=P(f)*H**hz[0]*Z**hz[1]
 ck('all original HZ coefficient equations',rec==den*v)
sc=load(str(p/'independent_subcritical_two_RR_delta_material_scaling.sobj'));U=sc['scale_ring'];rr,c,T,W,S=U.gens();UF=U.fraction_field();kh=K.ring().hom([c*c/rr,rr,UF.zero(),UF.zero()],UF);bh=sc['ring'].hom([c*c*rr,T,W,S],UF)
for it in sc['items']:
 f=it['original_source'];ck('coefficient source contains no base variables',all(co.derivative(H)==0 and co.derivative(Z)==0 for co in f.coefficients()));image=sum(kh(co.numerator())/kh(co.denominator())*(c*T)**int(ex[0])*W**int(ex[1])*(c*S)**int(ex[2]) for ex,co in f.dict().items());ck('all material scaling source identity',image==it['unit']*bh(it['primitive']))
curve=load(str(p/'independent_subcritical_two_RR_D_curve_sources.sobj'));Q=curve['ring'];N,X,y=Q.gens();BX=PolynomialRing(QQ,names=('check_N','check_T','check_y'));n,tt,yy=BX.gens();BXF=BX.fraction_field();h=sc['ring'].hom([n,tt,-yy/(1+yy),tt*yy],BXF);qh=Q.hom([n,tt**2,yy],BXF)
for it in curve['items']:
 pull=h(it['original_source']);saved=it['pullback'];savednew=BXF(BX({tuple(map(int,ex)):co for ex,co in saved.numerator().dict().items()}))/BXF(BX({tuple(map(int,ex)):co for ex,co in saved.denominator().dict().items()}));ck('exact original Dcurve pull',pull==savednew);re=it['rational_content']*tt**it['T_power']*yy**it['y_power']*(yy+1)**it['y_plus_one_power']*qh(it['primitive']);ck('complete source reconstruction after guarded divisions',pull.numerator()==re and all(z/(yy+1) in QQ for z,ep in pull.denominator().factor() if z not in QQ))
fact=load(str(p/'independent_subcritical_two_RR_D_curve_factor_ideals.sobj'));AA=fact['A_sources'];RR=fact['R_sources'];covered=set()
for src in curve['sources'][:-1]:
 match=[]
 for i,fa in enumerate(AA):
  for j,fr in enumerate(RR):
   rat=src/(fa*fr)
   if rat in QQ and rat:match.append((i,j));covered.add((i,j))
 ck('source is genuine product of A and R factors',bool(match))
ck('all eight products are necessary',covered==set((i,j) for i in range(4) for j in range(2)))
C=X*(y+1)**2+1;D=X*(y+1)**2*(y-1)+2*y+1;ck('original linear-even compatibility',fact['compatibility']/(C*D) in QQ and fact['compatibility'])
# Rebuild the original quadratic-even resultant itself from the original three e-equations.
evdata=load(str(p/'independent_subcritical_two_RR_even_resultants.sobj'));FL=evdata['flat_ring'];FF=FL.fraction_field();fa,fr,fH,fZ,fu,fw,ft=FL.gens();fh=K.ring().hom([fa,fr,fH,fZ],FF)
def flat(f):return sum(fh(co.numerator())/fh(co.denominator())*ft**int(ex[0])*fu**int(ex[1])*fw**int(ex[2]) for ex,co in f.dict().items())
quad=even[2];c2=quad.derivative(e,2)/2;c1=quad.derivative(e).subs({e:0});c0=quad.subs({e:0});ck('complete quadratic e source',quad==c2*e**2+c1*e+c0)
for i in range(2):
 A=even[i].derivative(e);B=even[i].subs({e:0});res=c2*B**2-c1*A*B+c0*A*A;ck('original necessary quadratic-even resultant '+str(i),flat(res)==evdata['items'][2]['items'][i]['flat_resultant'])
eq=load(str(p/'independent_subcritical_two_RR_D_curve_even_resultants.sobj'))
for it in eq['items']:
 prodval=FF(it['state_factor'])
 for z,pow in it['removed_factors']:prodval*=FF(z)**pow
 ck('all removed factors recorded',it['original_resultant']/prodval in QQ and it['original_resultant'])
 for z,pow in it['removed_factors']:ck('quadratic-even removed factor protected',any(FF(z)/FF(u) in QQ and FF(z)/FF(u) for u in [fa,fr,fa*fr+1,fH,fZ,fH+1,fZ-fH]))
cert=load(str(p/'independent_subcritical_two_RR_D_curve_exclusion_certificate.sobj'));ck('Rbranch full source guard identity',sum(c*f for c,f in zip(cert['R_cofactors'],cert['R_sources']))==cert['R_target']);ck('Rbranch protected target y square',cert['R_target']==y*y)
ap=load(str(p/'independent_subcritical_two_RR_D_curve_A_projection.sobj'));ak,ay=ap['ring'].gens();ah=ap['ring'].hom([N*X,y],Q)
for f,g in zip(ap['A_sources'],AA):ck('A branch source change K=NX',ah(f)==g)
ck('A branch exact univariate source lift',sum(c*f for c,f in zip(cert['A_cofactors'],cert['A_sources']))==cert['A_target']);ck('A projection polynomial exact',cert['A_target']==(ay**3-3*ay-1)*(ay**3+3*ay**2-1))
proj=load(str(p/'independent_subcritical_two_RR_D_curve_even_projection.sobj'));UP=proj['ring'];uy=UP.gen();for_names=[]
for it,bb in zip(proj['items'],cert['A_even_Bezout']):
 h=eq['ring'].hom([it['X_image'],uy],UP.fraction_field());vals=[h(z['primitive']) for z in eq['items']];ck('both original even pullbacks on compatibility chart',vals==it['values']);ck('gcd computed from both true even numerators',gcd(z.numerator() for z in vals)==it['gcd']);ck('removed y factor protected',bb['even_y']*uy**2==it['gcd']);ck('complete Bezout obstruction',bb['cofactors'][0]*bb['A_y']+bb['cofactors'][1]*bb['even_y']==1)
rec={'checks':len(checks),'seconds':time.monotonic()-st,'scope':'Representative two-RR subcritical Dcurve excluded, all complex directions and every material. Genuine free e is eliminated, never made base-constant.'};save({'checks':checks,'record':rec},str(p/'independent_subcritical_two_RR_D_curve_source_reader.sobj'),compress=False);(p/'independent_subcritical_two_RR_D_curve_source_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,flush=True)
