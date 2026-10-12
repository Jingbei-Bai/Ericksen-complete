"""Replay the T=2 all-spectral wall proof from the two full-Hb cubics."""
from pathlib import Path
import time,json
from sage.all import QQ,PolynomialRing,NumberField,load,save
out=Path(__file__).resolve().parent;start=time.monotonic();checks={}
src=load(str(out/'pencil_vzero_sone_cubics.sobj'))
B=PolynomialRing(QQ,names=('a','q'));a,q=B.gens();K=B.fraction_field()
R=PolynomialRing(K,'U');U=R.gen();r=q-1;L=a*q+1
AH=PolynomialRing(B,'Hb');Hb=AH.gen();full=[]
for p in src['cubics']:
 full.append([sum((B(c)*a**mon[0]*q**mon[1]*2**mon[2]*Hb**mon[3]
                  for mon,c in co.dict().items()),AH.zero()) for co in p.list()])
F=[]
for co in full[0]:
 quot,rem=co.quo_rem(Hb+1);assert rem==0;F.append(quot)
G=full[1]
checks['F_old_guard_content_only']=all(full[0][j]==(Hb+1)*F[j] for j in range(4))
series=[]
for cs,height in [(F,15),(G,20)]:
 series.append([R([co[height-i] for co in cs]) for i in range(4)])
F0,F1,F2,F3=series[0];G0,G1,G2,G3=series[1]
checks['F0_source']=F0==16*a*a*q*r**6*(U-L)*U*(2*r*U-r*L+4)
checks['G0_source']=G0==-16*a**3*q*q*r**7*U*(U-L)**2
checks['positive_spectral_identity']=1+r*L==q*(1+a*r)
fp=F0.derivative()(L);fs=-F1(L)/fp
gL1=G1(L);gL2=G2(L)+fs*G1.derivative()(L)+fs**2*G0.derivative(2)(L)/2
checks['L_simple_F_root']=fp==16*q*a*a*r**6*L*(r*L+4)
checks['L_first_zero']=gL1==0
checks['L_second_mismatch']=gL2==-16*q*a**3*r**9*L**3
gp=G0.derivative()(0);g1=-G1(0)/gp
g2=-(G2(0)+g1*G1.derivative()(0)+g1**2*G0.derivative(2)(0)/2)/gp
M1=F1(0)+g1*F0.derivative()(0)
M2=F2(0)+g1*F1.derivative()(0)+g2*F0.derivative()(0)+g1**2*F0.derivative(2)(0)/2
P=a*a*q*r**3+a*r*(r*r-12*r-8)-8*q
checks['zero_simple_G_root']=gp==-16*q*q*a**3*r**7*L**2
checks['zero_g1']=g1==-8/r
checks['zero_g2']=g2==24/r
checks['zero_first_mismatch']=M1==-128*q*a*r**3*L*P
checks['zero_second_is_polynomial']=M2.denominator()==1
CQ=PolynomialRing(QQ,'qq');qq=CQ.gen();AA=PolynomialRing(CQ.fraction_field(),'aa');aa=AA.gen()
def cvA(p):return sum((c*aa**mon[0]*qq**mon[1] for mon,c in p.dict().items()),AA.zero())
PP=cvA(P);MM=cvA(M2.numerator());RR=PP.resultant(MM)
checks['zero_resultant_exact']=RR==-8388608*qq**6*(qq-1)**12*(qq**2-4*qq+1)*(5*qq**2+508*qq+27)
remainder=MM.quo_rem(PP)[1]
records=[]
for index,qp in enumerate([qq**2-4*qq+1,5*qq**2+508*qq+27]):
 nf=NumberField(qp,'w');w=nf.gen()
 def evQ(x):return nf(x.numerator()(w))/nf(x.denominator()(w))
 cs=[evQ(c) for c in remainder.list()]
 checks[f'exception{index}_remainder_linear']=len(cs)==2 and cs[1]!=0
 av=-cs[0]/cs[1]
 checks[f'exception{index}_remainder_pivot_norm']=cs[1].norm()!=0
 if index==0:
  checks['exception0_a']=av==w-4
  checks['exception0_forces_L_zero']=av*w+1==0
 else:
  checks['exception1_a']=av==259*w/270+QQ(4871)/50
  def ev(x):return nf(x.numerator()(av,w))/nf(x.denominator()(av,w))
  NR=PolynomialRing(nf,'UU')
  ff=[NR([ev(c) for c in f.list()]) for f in series[0]]
  gg=[NR([ev(c) for c in g.list()]) for g in series[1]]
  f0,f1,f2,f3=ff;g0,g1p,g2p,g3p=gg
  den=g0.derivative()(0);checks['exception1_simple_root']=den!=0
  x1=-g1p(0)/den
  x2=-(g2p(0)+x1*g1p.derivative()(0)+x1*x1*g0.derivative(2)(0)/2)/den
  x3=-(g3p(0)+x1*g2p.derivative()(0)+x2*g1p.derivative()(0)
       +x1*x1*g1p.derivative(2)(0)/2+x1*x2*g0.derivative(2)(0)+x1**3*g0.derivative(3)(0)/6)/den
  m1=f1(0)+x1*f0.derivative()(0)
  m2=f2(0)+x1*f1.derivative()(0)+x2*f0.derivative()(0)+x1*x1*f0.derivative(2)(0)/2
  m3=f3(0)+x1*f2.derivative()(0)+x2*f1.derivative()(0)+x3*f0.derivative()(0)\
      +x1*x1*f1.derivative(2)(0)/2+x1*x2*f0.derivative(2)(0)+x1**3*f0.derivative(3)(0)/6
  checks['exception1_first_two_zero']=m1==m2==0
  checks['exception1_third_exact']=m3==-QQ(355260432384)/125*w-QQ(15897378816)/125
  checks['exception1_third_norm']=m3.norm()==QQ(107979287999067998650368)/15625
  records.append({'third_mismatch':str(m3),'third_norm':str(m3.norm())})
 records.append({'q_equation':str(qp),'remainder_pivot':str(cs[1]),
                 'remainder_pivot_norm':str(cs[1].norm()),'a':str(av)})
assert all(checks.values()),checks
payload={'source_file':'pencil_vzero_sone_cubics.sobj','series':series,
 'M1':M1,'M2':M2,'P':P,'resultant_a':RR,'remainder':remainder,'checks':checks,'records':records}
save(payload,str(out/'pencil_T2_uniform_audit.sobj'))
report={'status':'PASS_T2_ALL_POSITIVE_SPECTRA_WITH_aq_plus_1_NONZERO',
 'checks':checks,'count':len(checks),'records':records,'elapsed_seconds':time.monotonic()-start,
 'scope':'T=2; a>0, b=a*q/(1+a*(q-1))>0, simple spectrum; a*q+1 !=0. The crossing aq+1=0 is handled separately.'}
(out/'pencil_T2_uniform_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
