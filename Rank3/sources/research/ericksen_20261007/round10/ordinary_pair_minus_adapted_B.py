"""Generic-base exact adapted quadratic generators on m*l^2=-ar.
This is a polynomial source identity, not a normal rank inference.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
d=load(str(p/'ordinary_pair_double_kernel_coupled.sobj'));P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();q=1+r;k=a*r;m=1+k;G=H*Z-H+Z;l=P.gen();mu=m*l*l+k;M=d['matrix'];labs=d['labels']
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+650*1024**2,)*2)
def ck(n,x):assert bool(x),n;checks.append(n)
def rem(f):return sum((c*(-k/m)**(int(e)//2)*l**(int(e)%2) for e,c in P(f).dict().items()),P.zero())
def det(A):
    if not A.nrows():return P.one()
    i=min(range(A.nrows()),key=lambda j:sum(bool(v) for v in A.row(j)))
    return sum(((-1)**(i+j)*A[i,j]*det(A.delete_rows([i]).delete_columns([j])) for j in range(A.ncols()) if A[i,j]),P.zero())
rs=[labs.index(('old',i)) for i in (2,17,29)];cols=(1,2,6);J=M.matrix_from_rows_and_columns(rs,cols);D=det(J);Q=vector(P,[M[i,8] for i in rs]);N=[]
for j in range(3):
    JJ=matrix(P,J);JJ.set_column(j,Q);N.append(det(JJ))
ck('Cramer exact polynomial identity',J*vector(P,N)==D*Q)
assert D
out=[];first=[];failed2=[]
for i in range(M.nrows()):
    f=D*M[i,8]-sum((M[i,cols[j]]*N[j] for j in range(3)),P.zero());quo,res=f.quo_rem(mu**2)
    qq,rr=f.quo_rem(mu)
    ck('all original rows exact mu quotient '+str(labs[i]),not rr)
    out.append(quo);first.append(qq)
    if res:failed2.append(i)
DD=rem(D);nn=[rem(f) for f in N]
ck('adapted n1 at wall',nn[0]==m*DD)
ck('adapted n2 at wall',nn[1]==m*l*DD)
ck('adapted eta square at wall',nn[2]==0)
print('D wall degree',DD.degree(),'terms',len(DD.dict()),'rows not mu2',[labs[j] for j in failed2],flush=True)
for e,c in DD.dict().items():print('Dcoefficient',e,K(c).factor(),flush=True)
B_rows=[labs.index(('old',i)) for i in (1,2,17,29)]
B4=M.matrix_from_rows_and_columns(B_rows,(1,2,6,8));B=det(B4);bu=-384*H**2*Z**2/(r**3*a**3*(Z-H)**2*(Z-1)**2*m**2*G);bf=m*Z*(3*H+2)*l*l+k*(3*H+1)
ck('old full B determinant',B==bu*l*l*(m*l*l-k)*mu**2*bf)
newB=matrix(P,[[M[i,j] for j in cols]+[out[i]] for i in B_rows]);bn=det(newB)
ck('adapted B determinant exact',bn*mu**2==D*B)
ck('adapted B at wall exact protected form',rem(bn)==rem(D*bu*l*l*(m*l*l-k)*bf))
for i in failed2:
    print('NONZERO first quotient',labs[i],flush=True)
    for e,c in rem(first[i]).dict().items():print('coefficient',e,K(c).factor(),flush=True)
finalrows=[labs.index(('old',i)) for i in (2,17,29,35)];finalmatrix=[]
for i in finalrows:
    ck('B row whole u eta coefficient zero '+str(labs[i]),M[i,4]==0)
    ff=M[i,8]-m*M[i,1]-m*l*M[i,2];quo,rr=ff.quo_rem(mu)
    ck('direct first-shear mu quotient '+str(labs[i]),not rr)
    finalmatrix.append([M[i,j] for j in cols]+[quo])
finalmatrix=matrix(P,finalmatrix);finaldet=det(finalmatrix)
ck('new B4 wall only protected base factor',rem(finaldet)==-16*DD)
cyc=K.hom([a,r,K.one(),-K.one()],K)
ck('cyclic new B4 all materials',cyc(rem(finaldet)[0])==1536*k**2/m**3)
save({'ring':P,'source_matrix':M,'source_labels':labs,'pivot_rows':(2,17,29),'pivot_columns':cols,'pivot_matrix':J,'pivot_det':D,'adapted_numerators':N,'mu_quotients':first,'rows_not_mu_squared':failed2,'mu_squared_quotients':out,'B_rows':(1,2,17,29),'adapted_B_matrix':newB,'adapted_B_det':bn,'final_B_rows':(2,17,29,35),'final_B_matrix':finalmatrix,'final_B_det':finaldet,'checks':checks,'scope':__doc__},str(p/'ordinary_pair_minus_adapted_B.sobj'),compress=False)
rep={'status':'PASS WHOLE SOURCE QUADRATIC MODULE','checks':len(checks),'seconds':time.monotonic()-st}
(p/'ordinary_pair_minus_adapted_B.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep,flush=True)
