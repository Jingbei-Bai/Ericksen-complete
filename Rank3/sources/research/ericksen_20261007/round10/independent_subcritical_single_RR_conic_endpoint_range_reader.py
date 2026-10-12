"""Independent first/second-source replay for the conic endpoint closure.

Compare complete source polynomial derivatives with the old independent
original coefficient jets, check the actual six-source unit and norm,
and extract concrete odd/even second obstructions from original rows.
"""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_subcritical_single_RR_full_regular_chart_H2.sobj'));old=load(str(p/'independent_critical_single_RR3_normal.sobj'));qq=load(str(p/'independent_subcritical_single_RR_conic_endpoint_complete_quadratic.sobj'));T=qq['ring'];E=T.base_ring();K=E.base_ring();k,r=K.gens();xb=E.gen();v,s,w=T.gens();P=d['ring'];PK=P.base_ring();ph=PK.ring().hom([k,r],K);BK=old['ring'].base_ring();bh=BK.ring().hom([k/r,r,2,-1],K);base=(xb,-xb,E.zero(),E.zero(),E.zero(),2*xb,E.zero(),E.zero(),E.zero());checks=[]
def ck(n,val):assert bool(val),n;checks.append(n)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 n,dd=c.numerator(),c.denominator();gg=n.gcd(dd);n,dd=n//gg,dd//gg;lc=dd.leading_coefficient();return K(n/lc)/K(dd/lc)
def ne(c):return E([nk(u) for u in E(c).list()])
def em(f,hh,imgs):return ne(sum((hh(c.numerator())/hh(c.denominator())*prod(iv**int(n) for iv,n in zip(imgs,ee)) for ee,c in f.dict().items() if not any(n and not iv for n,iv in zip(ee,imgs))),E.zero()))
M=matrix(E,28,9,lambda i,j:em(d['sources'][i].derivative(P.gen(j)),ph,base))
for i in range(28):
 for j in range(9):ck('independent complete-source original normal '+str((i,j)),M[i,j]==em(old['normal_matrix'][i,j],bh,base[:7]))
R=matrix(E,9,6);R[0,0]=1;R[1,0]=-1;R[5,0]=2
for i,j in enumerate((3,4,6,7,8),1):R[j,i]=1
A=M*R;rows=(2,3,6,9,25,27);det=ne(A.matrix_from_rows(rows).determinant());ck('six actual source range minor nonzero',det!=0);print('DET6',det,flush=True)
mulmat=matrix(K,2,2,lambda i,j:ne(det*xb**j)[i]);norm=nk(mulmat.det());print('NORM6',norm.factor(),flush=True)
G=P.hom([P.gen(0),P.gen(1),P.zero(),P.zero(),P.zero(),P.gen(5),P.zero(),P.gen(7),P.zero()],P)
for i in range(25):ck('whole original odd source J5 identity '+str(i),G(d['sources'][i])==0)
left=qq['left_kernel'];Q=qq['quadratic_matrix'];obs=left*Q;gamma=E(-QQ(5)/504)*xb/(1+k);brow=next(i for i,row in enumerate(obs) if row[3] and not any(row[j] for j in (0,1,2,4,5)));scale=gamma/obs[brow,3];B=vector(E,[ne(scale*c) for c in left[brow]]);ck('odd residual uses only original odd sources',all(not B[i] for i in range(25,len(B))));ck('odd source range annihilation',vector(E,B[:28])*M==0);ck('literal odd second coefficient',vector(E,B)*Q==vector(E,(0,0,0,gamma,0,0)))
F=vector(E,[E.zero()]*len(qq['sources']));F[25]=1;F[26]=-1/k;ck('literal even range annihilation',vector(E,F[:28])*M==0);ck('literal even quadratic',F*Q==vector(E,(QQ(1)/18,QQ(7)/3,32,0,0,0)))
print('PASS',len(checks),'odd combination',[(i,c) for i,c in enumerate(B) if c],'sec',time.monotonic()-st,flush=True)
save({'ring':E,'normal_matrix':M,'range_columns':R,'range_rows':rows,'range_minor':det,'range_minor_norm':norm,'odd_combination':B,'odd_vw_unit':gamma,'even_combination':F,'checks':checks,'scope':__doc__},str(p/'independent_subcritical_single_RR_conic_endpoint_range_reader.sobj'),compress=False);(p/'independent_subcritical_single_RR_conic_endpoint_range_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st},indent=2)+'\n')
