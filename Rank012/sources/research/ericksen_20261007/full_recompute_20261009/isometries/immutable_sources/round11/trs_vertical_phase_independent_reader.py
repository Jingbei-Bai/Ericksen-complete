"""Independent replay of all-pitch finite-B -> finite-L bridge."""
from sage.all import *
from pathlib import Path
import json,time,resource
OUT=Path(__file__).resolve().parent;START=time.time();checks=[]
def ck(label,yes):
 assert yes,label
 checks.append(label)
first=load(str(OUT/'trs_algebra_regular_first_order_solution.sobj'));source=first['source']
K=source['field'];v,p,j,A,T,B,L=K.gens();sol=first['solution_coefficients'];M=source['M'];U2=source['U_squared']
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+650*1024**2,vm+650*1024**2))
for lab,column in sol.items():
 for row,z in enumerate(M*column+source['rhs'][lab]):ck('complete first-order source '+lab+' '+str(row),z==0)
dat=load(str(OUT/'trs_algebra_regular_frobenius.sobj'));bze,bzo=dat['B_zeta']
fac=2*B*L/(v+L*L)
def project(column):return 2*B*column[1]+fac*column[3]
expected_even={(0,0):project(sol['constant']),(1,0):project(sol['KC'])}
expected_odd={(0,0):project(sol['U'])+project(sol['one_over_U'])/U2,
              (1,0):project(sol['KC_over_U'])/U2,(0,1):project(sol['KB_over_U'])/U2}
for lab,actual,want in [('even',bze,expected_even),('odd',bzo,expected_odd)]:
 ck('complete '+lab+' pressure monomials',set(actual)==set(want))
 for mon,z in want.items():ck('direct source '+lab+' '+str(mon),actual[mon]==z)
ev=load(str(OUT/'trs_vertical_phase_even.sobj'));S=ev['ring'];SK=S.fraction_field();sv,sp,sj,sa,st,sb,sl,sc=map(SK,S.gens())
hom=K.hom([sv,sp,sj,sa,st,sb,sl],SK)
ck('full even numerator linked to original equations',SK(ev['numerator'])/ev['denominator']==hom(expected_even[(0,0)])+sc*hom(expected_even[(1,0)]))
n=ev['numerator'];ck('even numerator L degree at most three',n.degree(S.gen(6))<=3)
q0=n.coefficient({S.gen(6):0})-S.gen(0)*n.coefficient({S.gen(6):2})
q1=n.coefficient({S.gen(6):1})-S.gen(0)*n.coefficient({S.gen(6):3})
for idx,z in enumerate([q0,q1]):ck('pressure cancels before division '+str(idx),z.degree(S.gen(7))<=0)
old=load(str(OUT/'trs_vertical_phase_resultant.sobj'));F=old['field'];R=F.ring();vv,pp,jj,i1,i2,bb=map(F,R.gens())
P=vv**3-i1*vv*vv+i2*vv-jj;aa=jj/vv-P/(vv*bb);theta=vv*aa*(i2-aa-jj*(1+bb)/vv)-jj*vv*vv
polys=[]
for idx,poly in enumerate([q0,q1]):
 even=F(0);odd=F(0)
 for ex,co in poly.dict().items():
  z=co*vv**ex[0]*pp**ex[1]*jj**ex[2]*aa**ex[3]*theta**(ex[4]//2)*bb**ex[5]
  if ex[4]%2:odd+=z
  else:even+=z
 ck('exact required T parity '+str(idx),(even if idx==0 else odd)==0)
 z=(odd if idx==0 else even)/(aa**3*(vv-pp*pp)*(theta if idx else 1))
 ck('pitch cancels after protected division '+str(idx),z.derivative(pp)==0)
 nn=z.numerator();dd=z.denominator();gg=nn.gcd(dd);nn=nn//gg;dd=dd//gg
 ck('remaining denominator only powers v,B '+str(idx),len(dd.dict())==1 and all(all(ex[k]==0 for k in [1,2,3,4]) for ex in dd.dict()))
 # Numerical normalization is kept exactly as in source-field arithmetic.
 ck('reduced polynomial matches source archive '+str(idx),F(z.numerator())==F(old['polynomials'][idx]))
 polys.append(z.numerator())
res=polys[0].resultant(polys[1],R.gen(5));ck('entire resultant replay',res==old['resultant'])
Z=5*vv**3-3*i1*vv*vv+i2*vv+jj;H=F(res)/(524288*jj**6*vv**16*P**2*Z**3)
ck('H is polynomial',H.denominator()==1);H=R(H);ck('H full archive identity',H==old['H'])
SS=i2*i2-4*jj*i1
ck('H degree exactly eight',H.degree(R.gen(0))==8)
ck('H leading coefficient includes entire material wall',F(H.coefficient({R.gen(0):8}))==20*SS*SS)
ck('H next coefficient survives S=0',F(H.coefficient({R.gen(0):7}))==-132*i1*SS*SS-1804*jj*i2*SS-14112*jj**3)
report={'status':'PASS','checks':len(checks),'seconds':float(time.time()-START),
 'scope':'All pitch, all real pressures, all positive simple spectra; conditional on finite fixed-v B projection.',
 'protected_divisions':'A^3*(v-p^2)*T and A^3*(v-p^2)*T^2; remaining denominators are v and B powers.',
 'material_wall':'S=I2^2-4jI1=0 retained explicitly; H7=-14112j^3.',
 'labels':checks}
(OUT/'trs_vertical_phase_independent_reader.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='labels'},indent=2),flush=True)
