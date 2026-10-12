"""Independent full-source normal reconstruction and uniform algebraicity proof.
All spectral and horizontal-constant walls are covered by exact identities.
No Groebner basis, selected material value, or finite-fiber inference is used.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
old=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));pairs=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'))
K=pairs['coefficient_field'];a,r,H,Z=K.gens();q=1+r;k=a*r;m=1+k;G=H*Z-H+Z;N=a*r*r
O=old['ring'];OK=O.base_ring();ph=OK.hom([a,a*q/m,H/r,q*Z/r],K);cache={}
P=PolynomialRing(K,names=('l','u','v'));l,u,v=P.gens()
def ck(n,z):assert bool(z),n;checks.append(n)
def cv(c):
 c=OK(c)
 if c not in cache:cache[c]=ph(c)
 return cache[c]
def mp(f):return P({tuple(map(int,e)):cv(c) for e,c in O(f).dict().items()})
def source_row(lab):
 typ,i=lab
 if typ=='old':return [mp(old['matrix'][i,j]) for j in range(6)]
 it=next(z for z in pairs['items'] if z['pair']==i);f=it['full_homogeneous_pullback'];out=[]
 for j in range(6):
  ex=[0]*6;ex[j]=1
  out.append(P({tuple(map(int,e[:3])):K(c) for e,c in f.dict().items() if tuple(map(int,e[3:]))==tuple(ex)}))
 return out
def det(A):
 if not A.nrows():return P.one()
 i=min(range(A.nrows()),key=lambda i:sum(bool(z) for z in A.row(i)))
 return sum(((-1)**(i+j)*A[i,j]*det(A.delete_rows([i]).delete_columns([j])) for j in range(A.ncols()) if A[i,j]),P.zero())
sets=[(('old',1),('old',17),('old',5),('old',7),('old',11)),(('old',1),('old',17),('old',5),('old',7),('pair',2))]
units=[H/((H+1)**6*(Z-1)**4*r**4*a**3*(Z-H)*G**3*m**3),-H/(4*(H+1)**6*(Z-1)**4*r**4*a**3*(Z-H)*G**2*m**3)]
Fs=[];actual=[];matrices=[];r0=source_row(('old',0));ep=Z**2*G/(7*(Z-1)**2*(Z-H))
ck('actual eta0 protected unit times u',r0[4]==ep*u)
for j,rs in enumerate(sets):
 rows=[source_row(lab) for lab in rs];M=matrix(P,rows)
 ck('other eta columns exactly zero '+str(j),not any(M.column(4)))
 five=det(M.matrix_from_columns((0,1,2,3,5)));six=det(matrix(P,[r0]+rows))
 ck('true six determinant Laplace '+str(j),six==ep*u*five)
 F=P(five/(units[j]*(v+2)*l**(4-j)))
 ck('full source polynomial reconstruction '+str(j),five==units[j]*(v+2)*l**(4-j)*F)
 Fs.append(F);actual.append(six);matrices.append(matrix(P,[r0]+rows))
 print('SOURCE',j,'terms',len(F.dict()),'sec',time.monotonic()-st,flush=True)
R=PolynomialRing(K,names=('l','J','v'));ll,J,vv=R.gens();L=ll*ll;beta=H*(Z-1)/(Z*(H+1));subs=[]
for j,F in enumerate(Fs):
 f=R.zero()
 for e,c in F.dict().items():
  ck('even u parity '+str((j,e)),not e[1]%2)
  f+=c*ll**int(e[0])*(beta**2*J*(vv+2)*(vv-2*H))**(int(e[1])//2)*vv**int(e[2])
 subs.append(f)
def co(f,n):return R({(int(e[0]),int(e[1]),0):c for e,c in f.dict().items() if e[2]==n})
A0=-(7*H**2*Z+7*H**2+17*H*Z+6*H+6*Z+2)
B0=14*H**2*Z+14*H**2+13*H*Z+15*H+3*Z+4
R1=A0*(r*L+q)+B0*(J+q)
pref=k*m*m*H**5*(Z-1)**4/(H+1)*J*L*(m*L-k)*(m*L+k)
ck('generic leading v5 full identity',co(subs[0],5)==pref*R1)
ck('two base coefficient determinant',(-7)*3-14*(-6)==63)
W=PolynomialRing(K,names=('J','v'));wj,wv=W.gens()
def wall(f,lv,jv=None):
 ans=W.zero()
 for e,c in f.dict().items():
  ck('even l parity '+str(e),not e[0]%2)
  ans+=c*lv**(int(e[0])//2)*(wj if jv is None else jv)**int(e[1])*wv**int(e[2])
 return ans
ev=K.hom([a,r,K.one(),-K.one()],K)
def cyclic(f):return W({tuple(map(int,e)):ev(c) for e,c in f.dict().items()})
plus=cyclic(wall(subs[0],k/m))
def jc(f,dv,dj):return K(f.monomial_coefficient(wj**dj*wv**dv))
c0=[jc(plus,0,j) for j in (0,1)];c1=[jc(plus,1,j) for j in (0,1)];c2=[jc(plus,2,j) for j in (0,1)]
z1=c0[0]*c1[1]-c0[1]*c1[0];z2=c0[0]*c2[1]-c0[1]*c2[0]
f1=21622*N-989*(k+q);f2=14006*N+377*(k+q)
ck('plus wall first J-coefficient minor',z1==524288*k**8/m*f1)
ck('plus wall second J-coefficient minor',z2==262144*k**8/m*f2)
ck('plus wall protected Bezout',377*f1+989*f2==22003428*N)
minus=cyclic(wall(subs[0],-k/m))
ck('minus wall nonzero v constant coefficient',jc(minus,0,1)==-565248*k**4 and jc(minus,0,0)==0)
special=wall(subs[1],-q/r,-q);cv4=K(special.monomial_coefficient(wv**4));S=k+q
res=(-8*N-5*S)*H+(24*N+18*S)*(Z-1)
factor=-2/r**2*a*q**3*(Z-1)**4*H**4*(7*H*Z+7*H+3*Z+4)*S*m**2
ck('remaining wall v4 exact identity',cv4==factor*res)
ck('remaining wall base-coefficient Bezout',18*(-8*N-5*S)+5*(24*N+18*S)==-24*N)
ck('remaining wall S0 already minus relation',K((m*(-q/r)+k).subs({a:-(1+r)/r}))==0)
save({'ring':P,'source_matrices':matrices,'six_minors':actual,'factors':Fs,'factor_units':units,'first_integral_ring':R,'first_integral_polynomials':subs,'generic_v5_remainder':R1,'plus_wall_coefficient_minors':(z1,z2),'plus_wall_primitive_factors':(f1,f2),'minus_wall_constant':-565248*k**4,'special_wall_v4':cv4,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_full_algebraicity_reader.sobj'),compress=False)
rec={'status':'PASS ALL MATERIAL ALL NONZERO L AND J ALGEBRAICITY','checks':len(checks),'seconds':time.monotonic()-st,'remaining':'This proves v algebraic off v=-2, v=2H and u=0; u=0 separately certified. The two sheets and l=0 require no-lift proofs.'}
(p/'ordinary_pair_full_algebraicity_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print(rec,flush=True)
