"""Replay the full positive-spectrum, real-pressure L=infinity certificate.

No sample spectra, floating point arithmetic, or material-dependent division.
"""
from sage.all import *
from pathlib import Path
import json,time
OUT=Path(__file__).resolve().parent; START=time.time();checks=[]
def ck(label,condition):
 if not condition:raise AssertionError(label)
 checks.append(label)

# Exact source identity. Bcompat has no new equation or new source hypothesis.
su=load(str(OUT/'tr_pure_rotation_logU_frobenius.sobj'))
sl=load(str(OUT/'tr_pure_rotation_L_frobenius.sobj'))
sb=load(str(OUT/'tr_pure_rotation_Bcompat_frobenius.sobj'))
S=su['P'];KS=S.fraction_field();sv,sj,sA,sT,sB,sL=map(KS,S.gens())
for lab in su['coefficients']:
 ck('source Bcompat '+lab,KS(sb['coefficients'][lab])==KS(sl['coefficients'][lab])+(sv+sL*sL)/sL*KS(su['coefficients'][lab]))

# Recompute every coefficient in the finite Newton support from the archives.
src=OUT/'tr_pure_rotation_B_infinity_leaders.py'
ns={'__file__':str(src),'__name__':'leader_replay'}
exec(src.read_text().split('vm=int(')[0],ns)
lead=ns['leader'];LL=ns['L']
sp=load(str(OUT/'tr_infinite_angle_newton_support.sobj'))
for name,dat in [('L_frobenius',sl),('Bcompat_frobenius',sb)]:
 for lab,z in dat['coefficients'].items():
  for which,poly in [('numerator',z.numerator()),('denominator',z.denominator())]:
   prof=sp['profiles'][name][lab][which];by={}
   for mon,co in poly.dict().items():by.setdefault(mon[5],{})[mon]=co
   ck(name+' '+lab+' '+which+' support',set(prof)==set(by))
   for ell,dic in by.items():
    n,c=lead(poly.parent()(dic),60)
    ck(name+' '+lab+' '+which+' coefficient '+str(ell),(n,c/LL**ell)==prof[ell])
   # Linear inequalities certify both complete open exponent intervals.
   for interval,r in [('0<r<1',QQ(1)/2),('r>1',QQ(2))]:
    low=min(n-r*b for b,(n,c) in prof.items())
    wins=[b for b,(n,c) in prof.items() if n-r*b==low]
    ck(name+' '+lab+' '+which+' unique cone leader '+interval,len(wins)==1)
    win=wins[0];nwin=prof[win][0]
    for b,(n,c) in prof.items():
     a=n-nwin; slope=win-b
     ck(name+' '+lab+' '+which+' full cone '+interval+' '+str(b),
        (a>=0 and a+slope>=0) if r<1 else (a+slope>=0 and slope>=0))
   if which=='denominator':
    # At r=1 the denominator also has one surviving L monomial.
    lo=min(n-b for b,(n,c) in prof.items());ws=[b for b,(n,c) in prof.items() if n-b==lo]
    ck(name+' '+lab+' denominator at r=1 single monomial',len(ws)==1)
    co=prof[ws[0]][1];nn=co.numerator();dd=co.denominator()
    ck(name+' '+lab+' denominator protected',len(nn.dict())==len(dd.dict())==1 and
       all(tuple(mon[2:])==(0,0,0) for pol in [nn,dd] for mon in pol.dict()))

# The balanced-angle resultant and its nonvanishing top v coefficient.
pd=load(str(OUT/'tr_infinite_angle_primitive_polynomials.sobj'))
f,g=pd['f'],pd['g'];K0=f.base_ring();P=K0.ring();v,j,I1,I2,C,D=P.gens();ii=P.base_ring().gen()
Kprof=sp['K'];RR=PolynomialRing(P.base_ring(),names=('w','j','I1','I2','C','D','ell'));KK=RR.fraction_field()
ww,jj,i11,i22,cc,dd,ell=map(KK,RR.gens());vv=ww*ww
hp=Kprof.hom([ww,jj,vv**3-i11*vv**2+i22*vv-jj,i22,ell],KK)
hs=K0.hom([vv,jj,i11,i22,cc,dd],KK)
def leading_sum(name,r):
 terms={}
 for lab,record in sp['profiles'][name].items():
  nd=[]
  for which in ['numerator','denominator']:
   prof=record[which];order=min(n-r*b for b,(n,c) in prof.items())
   coeff=sum(hp(c)*ell**b for b,(n,c) in prof.items() if n-r*b==order);nd.append((order,coeff))
  order=nd[0][0]-nd[1][0]+(1-r if lab in ['EC','ED'] else 0)
  factor={'1':1,'EC':ww*ell*cc,'ED':ww*ell*dd,'CC':cc*cc,'CD':cc*dd,'DD':dd*dd}[lab]
  terms[order]=terms.get(order,KK(0))+nd[0][1]/nd[1][1]*factor
 lo=min(terms);return lo,terms[lo]
for name,poly,order in [('L_frobenius',f,2),('Bcompat_frobenius',g,4)]:
 lo,z=leading_sum(name,QQ(1));expect=sum(hs(co)*(ell/(ii*ww))**k for k,co in enumerate(poly.list()))
 ratio1=z/expect;nn=ratio1.numerator();de=ratio1.denominator();gg=nn.gcd(de);nn=nn//gg;de=de//gg
 ck(name+' balanced exponent',lo==order)
 ck(name+' balanced source agrees polynomial up to protected monomial',len(nn.dict())==len(de.dict())==1 and
    all(all(mon[k]==0 for k in [2,3,4,5]) for pp in [nn,de] for mon in pp.dict()))
lo,z=leading_sum('L_frobenius',QQ(1)/2)
ck('all 0<r<1 leading expression',lo==2 and z==vv/jj*(6*i11-7*ii*cc))
H2=vv**2*(i11**2-QQ(5)/2*ii*i11*cc-cc**2-ii*jj*dd+2*i22)+vv*(-ii*jj*i11*dd-jj*cc*dd-6*jj)+jj*(-3*i11+QQ(7)/2*ii*cc)
lo,z=leading_sum('L_frobenius',QQ(2))
ck('all r>1 leading expression',lo==0 and z==2*ell**2/jj**2*H2)
raw=load(str(OUT/'tr_infinite_angle_leading_polynomials.sobj'))
ck('primitive f only constant scaling',f==raw['f']*pd['f_scale'])
ck('primitive g only constant scaling',g==raw['g']*pd['g_scale'])
res=f.resultant(g)
ck('resultant polynomial',res.denominator().total_degree()==0)
res=P(res);ck('resultant v degree six',res.degree(v)==6)
ck('resultant divisible v fourth',all(res.coefficient({v:k})==0 for k in range(4)))
F=load(str(OUT/'tr_infinite_angle_top_resultant_real_0.sobj'))['f']
G=load(str(OUT/'tr_infinite_angle_top_resultant_real_1.sobj'))['f']
ratio=K0(res.coefficient({v:6}))/(F*G)
ck('highest resultant coefficient F G protected',ratio==K0(ratio.numerator())/ratio.denominator() and
   len(ratio.numerator().dict())==len(ratio.denominator().dict())==1 and
   all(mon[0]==mon[2]==mon[3]==mon[4]==mon[5]==0 for poly in [ratio.numerator(),ratio.denominator()] for mon in poly.dict()))
def realpart(p,k):return P({m:QQ(co[k]) for m,co in p.dict().items()})
Fr,Fi=realpart(F,0),realpart(F,1)
positive_bracket=I1*(C+3*j*D/I1)**2+j*(I2-9*j/I1)*D**2+12*(I1*I2-9*j)
ck('F imaginary positive square identity',K0(Fi)==4*C*positive_bracket)
ck('F real at C zero strict sign',Fr.subs({C:0})==-36*(j*j*D*D+4*(I2*I2-3*j*I1)))
Gr,Gi=realpart(G,0),realpart(G,1);rd=Gr.resultant(Gi,D)
ck('G real resultant agrees archive',rd==load(str(OUT/'tr_infinite_angle_top_resultant_real_1.sobj'))['resultant'])

# Normalize by I1>0: x=I2/I1² in (0,1/3], y=j/I1³ in (0,1/27].
R=PolynomialRing(QQ,names=('x','y','z'));x,y,z=R.gens();q=R(0)
rr=P(rd)/(64*j*j)
ck('G real resultant normalization polynomial',rr.denominator().total_degree()==0)
rr=P(rr)
for mon,co in rr.dict().items():
 vv,jj,i1,i2,cc,dd=mon
 ck('G resultant even C and homogeneous '+str(mon),vv==dd==0 and cc%2==0 and 3*jj+i1+2*i2+cc==14)
 q+=QQ(co)*x**i2*y**jj*z**(cc//2)
ck('G positive constant coefficient',q.coefficient({z:0})==y+2*x*y+QQ(9)/4*y*y)
bounds={}
for k in range(1,8):
 pol=q.coefficient({z:k});bound=QQ(0)
 for mon,co in pol.dict().items():
  a,b,zz=mon
  if a==b==0:bound+=co
  elif co<0:bound+=co*(QQ(1)/3)**a*(QQ(1)/27)**b
 ck('G strictly positive z coefficient '+str(k),bound>0)
 bounds[str(k)]=str(bound)

report={'status':'PASS','checks':len(checks),'seconds':float(round(time.time()-START,3)),
 'scope':'All positive simple spectra, all real KC and KB, every Puiseux growth L~ell tau^-r with r>0; outside finitely many v values.',
 'resultant_top_ratio':str(ratio),'positive_coefficient_lower_bounds':bounds,
 'exceptions':'Only finitely many v levels; no pressure specialization or omitted growth exponent.',
 'checks_labels':checks}
(OUT/'tr_infinite_angle_whole_reader.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
save({'resultant':res,'F':F,'G':G,'G_normalized_real_resultant':q,'bounds':bounds,'profile_archive':'tr_infinite_angle_newton_support.sobj'},str(OUT/'tr_infinite_angle_whole_certificate.sobj'))
print(json.dumps({k:v for k,v in report.items() if k!='checks_labels'},indent=2),flush=True)
