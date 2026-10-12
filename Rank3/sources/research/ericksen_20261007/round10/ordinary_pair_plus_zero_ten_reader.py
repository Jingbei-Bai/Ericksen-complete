"""Exact whole-source ten-generator certificate on u=0,v=-2,m*l^2=ar.
The source ideal is relative in u and v; neither is frozen in the
membership statement.  In particular the pure t^2 divisibility is
checked at arbitrary u before evaluating the coefficient matrix.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(n,x):assert bool(x),n;checks.append(n)
d=load(str(p/'ordinary_pair_double_kernel_coupled.sobj'));cr=load(str(p/'ordinary_pair_plus_wall_cubic_reader.sobj'));q2=load(str(p/'independent_incoming_simple_wall_quadratics.sobj'));nr=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));ep=load(str(p/'ordinary_pair_plus_zero_E29_cubic.sobj'))
P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();q=1+r;k=a*r;m=1+k;G=H*Z-H+Z;l=P.gen();mu=m*l*l-k;L=l*l
OK=q2['ring'].base_ring();ph=OK.hom([a,a*q/m,H/r,q*Z/r],K);OR=q2['ring'];ol,ou,ov=OR.gens();U2=PolynomialRing(K,names=('l','u'));ll,uu=U2.gens()
def fullmap(f):return sum((ph(c)*ll**int(e[0])*uu**int(e[1])*(-2)**int(e[2]) for e,c in OR(f).dict().items()),U2.zero())
pure=q2['items'][1]['values'];assert q2['items'][1]['name']=='k'
bs=[('old',i) for i in (1,2,17,29)];ars=[('old',i) for i in (5,7,11,13)]+[('pair',2)];labs=[('old',0)]+bs+ars
for lab in labs:
    if lab[0]!='old':continue
    i=lab[1];pred=(8*(m*ll*ll-k) if i==2 else -48*ll*(m*ll*ll-k) if i==29 else U2.zero())
    ck('all-u pure quadratic from complete twojet '+str(i),fullmap(pure[i]['second'])==pred)
    ck('whole eta coefficient divisible u '+str(i),not nr['matrix'][i,4].subs({ou:0}))
    ck('whole kernel coefficient divisible v+2 '+str(i),not (nr['matrix'][i,3]-nr['matrix'][i,0]).subs({ov:-2}))
pairs=load(str(p/'ordinary_pair_incoming_full_twojet.sobj'));pair2=pairs['items'][2];PR=pairs['jet_ring'];pl,pu,pv=PR.gens()
ck('pair2 whole eta coefficient arbitrary v',pair2['jet_coefficients'][5]==0)
ck('pair2 whole kernel coefficient arbitrary u',not (pair2['jet_coefficients'][4]-pair2['jet_coefficients'][1]).subs({pv:-2}))
# Independently extract H7 (=old17) from its full signed products on the
# pure kernel ray, keeping the ordinary total degree through the product.
hs=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));RO=PolynomialRing(OK,'l');rl=RO.gen();cache={}
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
def factor(i):
    if i in cache:return cache[i]
    ans={}
    for exp,c in hs['polynomial_pool'][i].items():
        e=tuple(map(int,exp));j=e[1]+e[4];degree=sum(e)
        if e[2] or e[3] or e[7] or j>3 or degree>5:continue
        ans[(degree,j)]=ans.get((degree,j),RO.zero())+c*(-1)**e[1]*(-2)**(e[6]+e[8])*rl**(e[5]+e[6])
    cache[i]=ans;return ans
h17=RO.zero()
for sign,ia,ib in hs['records'][7]['full_signed_product_DAG']:
    A=factor(ia);B=factor(ib)
    h17+=sign*sum((f*B.get((5-deg,3-j),RO.zero()) for (deg,j),f in A.items()),RO.zero())
ck('full signed H7 cubic zero',not h17)
ck('full original E29 cubic zero',not ep['full_pure_coefficients'][3])
ec2=sum((ph(c)*l**int(e) for e,c in ep['full_pure_coefficients'][2].dict().items()),P.zero())
ck('full E29 direct source quadratic',ec2==-48*l*mu)
def cub(lab):
    if lab in (('old',17),('old',29)):return P.zero()
    ff=cr['full_pure_kernel_restrictions'][lab][3]
    return sum((K(c)*l**int(e[0]) for e,c in ff.dict().items() if not e[1]),P.zero())
rows=[]
for lab in labs:
    row=list(d['matrix'].row(d['labels'].index(lab)));quo,remainder=row[8].quo_rem(mu)
    ck('whole polynomial square quotient '+str(lab),not remainder)
    rows.append(row[:8]+[quo,cub(lab)])
M=matrix(P,rows)
ac=(0,3,5,7,9);bc=(1,2,6,8);perm=(4,)+bc+ac
for j in range(1,5):ck('complete B support '+str(j),not any(M[j,c] for c in (4,)+ac))
for j in range(5,10):ck('complete A support '+str(j),not any(M[j,c] for c in (4,)+bc))
U=Z**2*G/(7*(Z-1)**2*(Z-H));ck('only row0 ueta unit',M[0,4]==U)
def det(A):
    if not A.nrows():return P.one()
    i=min(range(A.nrows()),key=lambda j:sum(bool(v) for v in A.row(j)))
    return sum(((-1)**(i+j)*A[i,j]*det(A.delete_rows([i]).delete_columns([j])) for j in range(A.ncols()) if A[i,j]),P.zero())
A=det(M.matrix_from_rows_and_columns(range(5,10),ac));B=det(M.matrix_from_rows_and_columns(range(1,5),bc));sgn=Permutation([j+1 for j in perm]).signature()
ck('full true ten determinant',det(M)==sgn*U*A*B)
bu=-384*H**2*Z**2/(r**3*a**3*(Z-H)**2*(Z-1)**2*m**2*G)
bf=m*Z*(3*H+2)*L+k*(3*H+1)
ck('whole B polynomial factorization',B==bu*L*(m*L+k)**2*bf)
def rem(f):return sum((c*(k/m)**(int(e)//2)*l**(int(e)%2) for e,c in P(f).dict().items()),P.zero())
D=k*r+m*q*Z**2
F=-11*H*Z**2+Z**3-12*H*Z-4*Z**2-H+2*Z+1
au=32*a**4*r**3*H**2*Z*(Z-1)/((Z-H)**2*(H+1)**2*m**4*G**2)
ck('A wall exact material/base factorization',rem(A)==au*D*F)
ck('B wall exact material/base factorization',rem(B)==bu*(k/m)*(2*k)**2*k*(Z*(3*H+2)+3*H+1))
ck('D generic base leading unit',D==k*r+m*q*Z**2)
ck('F nonzero base polynomial',K(F).derivative(H)==-(11*Z**2+12*Z+1))
ck('last B base factor nonzero',K(Z*(3*H+2)+3*H+1).derivative(H)==3*(Z+1))
save({'ring':P,'labels':labs,'matrix':M,'A_rows':ars,'A_columns':ac,'A_determinant':A,'B_rows':bs,'B_columns':bc,'B_determinant':B,'A_wall_unit':au,'A_wall_base_factors':(D,F),'B_unit':bu,'B_base_factor':bf,'column_permutation':perm,'permutation_sign':sgn,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_plus_zero_ten_reader.sobj'),compress=False)
report={'status':'PASS WHOLE POSITIVE AMPLITUDE WALL AT u=0','checks':len(checks),'seconds':time.monotonic()-st,'scope':'True complete relative J10 ideal, all u/v corrections retained; all materials a*r*q*m!=0, generic base, arbitrary ramification.'}
(p/'ordinary_pair_plus_zero_ten_reader.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)
