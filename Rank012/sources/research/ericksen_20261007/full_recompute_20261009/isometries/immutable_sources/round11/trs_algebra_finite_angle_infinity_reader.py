"""Source audit of finite nonzero L at B=infinity and its two exceptional centers."""
exec(open(__file__.replace('trs_algebra_finite_angle_infinity_reader.py','trs_algebra_frobenius_infinity_leading.py')).read().split('items=[]')[0])
checks=[]
def ck(name,cond):
 if not cond:raise AssertionError(name)
 checks.append(name)
lead=load(str(OUT/'trs_algebra_frobenius_infinity_leading.sobj'))
rlead=lead['base'].ring();mapping=lead['base'].hom([v,p,j,I2,P,L],F)
# Coefficients of 1,Y, where Y=(v/j)*tau*u and Y²=-v/(Q D).
lc={}
for item in lead['items']:
 if item['twice_B_power']!=-2:continue
 odd=item['U_parity'];mon=tuple(item['pressure_monomial']);row=item['row']
 ck('root_parity '+str((row,odd,mon)),item['tau_parity']==odd)
 coef=mapping(item['coefficient'])*(j/v if odd else 1)
 lc[row,odd,mon]=nf(coef)
q=v-p*p;D=L*L+v;M=p*p*L*L+2*v*p*p-v*v
I1f=(v**3+v*I2-j-P)/v**2
bE=-p*(p*p*L*L+4*v*p*p-3*v*v)
bO=-q*(3*p*p*L*L+4*v*p*p-v*v)
for odd,target in [(0,bE),(1,bO)]:
 ck('CB_constant '+str(odd),nf(lc[0,odd,(0,0)]+L/D*lc[1,odd,(0,0)])==0)
 ck('CB_KC '+str(odd),nf(lc[0,odd,(1,0)]+L/D*lc[1,odd,(1,0)]-target/(j*L*L*q))==0)
ck('C1_constant_even',nf(lc[1,0,(0,0)]+6*I1f*M/(j*v))==0)
ck('C1_constant_odd',nf(lc[1,1,(0,0)]+12*I1f*p*q*D/(j*v))==0)
ck('KC_nonzero_norm',nf(q*D*bE*bE+v*bO*bO-q*(p*p*L*L+v*v)**3)==0)
ck('KC_zero_norm',nf(M*M+4*p*p*q*v*D-(p*p*L*L+v*v)**2)==0)
# At any finite L with L!=0,D!=0, all source denominators are z-normalized units.
allowed=[v,j,L,p*p-v,v-p*p,D]
for row,C in enumerate(d['frobenius']):
 for odd,pol in enumerate(C):
  for mon,x in sorted(pol.items()):
   nn=poly_lead(x.numerator());dd=poly_lead(x.denominator())
   ck('z_regular '+str((row,odd,mon)),-nn[0]+dd[0]-odd>=2)
   ck('protected_finite_L_denominator '+str((row,odd,mon)),all(nf(ff) in allowed for ff,ee in dd[1].factor()))
# Ideal membership in the analytic local ring, using the formal free h=delta/z.
# A homogeneous degree-m Taylor polynomial becomes z^m times a polynomial in h;
# hence identically zero h-polynomials at degrees 0,1 check every such monomial.
blo=load(str(OUT/'trs_algebra_newton_blowup_1.sobj'));BK=blo['base']
bv,bp,bj,bI1,bI2,bC,bD,bh=BK.gens();bii=BK.ring().base_ring().gen();bhp=BK.ring().gen(7)
cs=blo['residuals']['C1']
ck('analytic_degree_zero_vanishes',cs[2]==0)
ck('analytic_degree_one_entire_polynomial_vanishes',cs[3]==0)
c4=BK(cs[4]);ck('degree_two_face_polynomial',c4.numerator().degree(bhp)==2 and c4.denominator().degree(bhp)==0)
kappa=BK(c4.numerator().coefficient({bhp:2}))/c4.denominator()
ck('delta_squared_protected_unit',kappa==-bv*bp*bp*(6*bI1+7*bii*bC)/(bj**3*(bp*bp-bv)))

# Strict second-face source and positivity evidence are independently replayed by these readers.
face=load(str(OUT/'trs_algebra_newton_face_polynomials.sobj'))
FK=face['base'];fv,fp,fj,fI1,fI2,fC,fD,fh=FK.gens()
Qf,Rf=face['face_polynomials']
ck('Q_h_degree',Qf.degree(FK.ring().gen(7))==2)
ck('Q_v_free',Qf.degree(FK.ring().gen(0))==0)
ck('R_v_linear',Rf.degree(FK.ring().gen(0))==1)
ck('Q_protected_h_leading',FK(Qf.coefficient({FK.ring().gen(7):2}))==-fj*(fI1+FK.ring().base_ring().gen()*QQ(7)/6*fC))
report={'status':'PASS_FINITE_NONZERO_L_INFINITY_SOURCE','checks':len(checks),'seconds':time.time()-START,
'centers':'L=+/- i v/p; both signs related by complex conjugation',
'newton_source_reader':'trs_algebra_newton_full_source_reader.py: 13 exact full-source checks',
'nonzero_resultant_for_every_real_parameters':'trs_newton_face_positive_transfer_reader.py transfers the 507-check all-parameter F/G positivity theorem',
'scope':'Generic v for each fixed positive simple spectrum and arbitrary real p!=0,KC,KB. Does not handle L=0,D=0,L=infinity.'}
(OUT/'trs_algebra_finite_angle_infinity_reader.json').write_text(json.dumps(report,indent=2));signal.alarm(0);print(json.dumps(report),flush=True)

