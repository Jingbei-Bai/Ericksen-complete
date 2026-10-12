"""Independent five-normal replay from the complete missing-RR1 chart.

The old source identity is polynomial in t, hence includes t=0. Check
every 25x5 transported coefficient, then an independent QQ Bezout unit.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_three_p_five_module.sobj'));E=d['ring'];K=E.base_ring();k,r=K.gens();m=1+k;q=1+r;x=E.gen();u=x*x/m;z=x*(1-u);olde=load(str(p/'independent_subcritical_single_RR_three_p_second_normal.sobj'))['candidate_images'][5];old=load(str(p/'independent_subcritical_two_RR_full_normal_case0.sobj'));OP=old['ring'];OK=OP.base_ring();ph=OK.ring().hom([k/r,r,2,-1],K);checks=[];cache={}
def ck(n,v):assert bool(v),n;checks.append(n)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 n,den=c.numerator(),c.denominator();gg=n.gcd(den);n,den=n//gg,den//gg;lc=den.leading_coefficient();return K(n/lc)/K(den/lc)
def ne(v):return E([nk(v[i]) for i in range(8)])
Q=PolynomialRing(QQ,'fr_u');v=Q.gen();F=v**4-2*v**3+2*v*v-v+4;P=12*v**3+123*v*v-473*v+468;gcd,invlin,unused=(1-v).xgcd(F);assert gcd==1;h=Q.hom([u],E);xinv=-(x**7-2*m*x**5+2*m*m*x**3-m**3*x)/(4*m**4);zinv=ne(xinv*h(invlin));ck('protected candidate z inverse',ne(z*zinv)==1);images=(E.zero(),ne(x*zinv),zinv,ne(olde*zinv*zinv))
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=nk(ph(c.numerator())/ph(c.denominator()))
 return cache[c]
def mp(f):
 ans=E.zero()
 for ex,c in f.dict().items():
  ex=tuple(map(int,ex))
  if ex[0]:continue
  term=E(cv(c))
  for image,n in zip(images,ex):
   for unused in range(n):term=ne(term*image)
  ans=ne(ans+term)
 return ans
N=d['matrix'];transport=matrix(E,25,5)
for i in range(25):
 degree=2 if i==6 else 4;M=old['matrix'];transport[i,0]=ne(z**degree*mp(M[i,4]));transport[i,1]=ne(z**(degree-2)*mp(M[i,0]));transport[i,2]=ne(z**(degree-2)*mp(M[i,1]));transport[i,3]=ne(z**(degree-3)*mp(M[i,3]));transport[i,4]=ne(z**(degree+1)*mp(M[i,2]))
 for j in range(5):ck('independent full source normal '+str((i,j)),transport[i,j]==N[i,j])
for i in range(25):
 label=('R boundary zero '+str(i)) if i<6 else ('Q0 boundary zero' if i==6 else 'H complete boundary zero '+str(i-7));ck('ambient whole-t source zero identity '+str(i),label in old['checks'])
cert=load(str(p/'independent_subcritical_single_RR_three_p_five_minor_2_3_6_21_9.sobj'));det=cert['determinant'];ck('exact candidate determinant polynomial',det==E(QQ(16)/42525*k*r*r*q*m**3)*h(P));gg,bF,bP=F.xgcd(P);ck('independent QQ polynomial unit',gg==1 and bF*F+bP*P==1);ck('all eight complex roots covered',F.discriminant()!=0 and F(0)!=0 and F(1)!=0 and F(2)!=0)
print('PASS',len(checks),'sec',time.monotonic()-st,flush=True);save({'ring':E,'transported_matrix':transport,'QQ_Bezout':(bF,bP),'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_three_p_five_reader.sobj'),compress=False);(p/'independent_subcritical_single_RR_three_p_five_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n')
