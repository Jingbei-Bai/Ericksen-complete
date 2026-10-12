"""Root replay: original complete Q6, original h/z derivatives, four-root quotient.

Build the complete differential numerator first in the original material
field, then specialize its coefficients. No array-jet source routine or
transformed-Q coefficient derivative is reused from the creator.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
def ck(s,b):assert b,s;checks.append(s)
g=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'));raw=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));wall=load(str(p/'independent_subcritical_single_RR_conic_wall.sobj'));ref=load(str(p/'independent_subcritical_single_RR_complete_E_wall_reader.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+750*1024**2,)*2)
O=g['ring'];OK=O.base_ring();oa,ob,oh,oz=OK.gens();ov=O.gens();source=O({tuple(map(int,e)):OK(c) for e,c in raw['quartic_source_rows'][6].dict().items()});it=g['items'][1]
parts=[(f,source.derivative(ov[j])) for j,f in enumerate(it['field_numerators'])]
for j,bp in enumerate((oh,oz)):parts.append((it['base_numerators'][j],O({tuple(map(int,e)):c.derivative(bp) for e,c in source.dict().items()})))
K=PolynomialRing(QQ,names=('rootwall_k','rootwall_r')).fraction_field();k,r=K.gens();m=1+k;q=1+r;ph=OK.ring().hom([k/r,k*q/(r*m),2/r,-q/r],K);cc={}
def cv(co):
 co=OK(co)
 if co not in cc:cc[co]=ph(co.numerator())/ph(co.denominator())
 return cc[co]
R=PolynomialRing(K,names=O.variable_names());rv=R.gens()
def cp(f):return R({tuple(map(int,e)):cv(co) for e,co in f.dict().items()})
G=sum((cp(f)*cp(df) for f,df in parts),R(0));ck('complete G ordinary degree bound',G.total_degree()<=6)
danger=R({e:co for e,co in G.dict().items() if sum(e)==6 and sum(e[6:])==3});ck('complete weighted9 face exactly zero',danger==0)
print('COMPLETE G TERMS',len(G.dict()),'seconds',time.monotonic()-st,flush=True)
A=PolynomialRing(K,names=('rootW','rootX','rootKap','rootT'));aw,ax,ak,at=A.gens();J=A.ideal(aw**2-k,ax**2-4*m,ak**3,at**2);B=A.quotient(J,names=('W','X','Kap','T'));W,X,Kap,T=B.gens()
EW=wall['ring'];wk=EW.base_ring().base_ring();wh=wk.ring().hom([k,r],K)
def wc(z):return wh(z.numerator())/wh(z.denominator())
def root(z):
 z=EW(z);return B(sum((wc(z[i][j])*X**i*W**j for i in range(2) for j in range(2)),B(0)))
ker=[root(z) for z in wall['kernel'][0]];ck('formal C normal kernel zero',ker[4]==0);ck('formal s2 normal kernel zero',ker[5]==0);ck('actual free t and no kappa',ker[6]==1 and ker[7]==0)
e0=-2*X*W+ker[0]*T;e1=-W/2+ker[1]*T;e2=2*X+ker[2]*T;al=B(1);bb=B(QQ(4)/3);de=ker[3]*T
images=[X,B(1),-X,W,-X*T/bb-Kap*e2/bb,X*W*T/bb+Kap*(X*e1-e0)/al,al*W*T/bb+Kap*e1,T,B(1)];hp=R.hom(images,B)
ck('gauge denominator chart unit',hp(cp(it['denominator']))==1)
ans={}
for dg,rd in sorted(set((sum(e),sum(e[6:])) for e in G.dict())):
 dp=6-dg;kp=13-2*dg-rd;ck('nonnegative delta weight '+str((dg,rd)),dp>=0)
 if dp>=2:continue
 f=R({e:co for e,co in G.dict().items() if sum(e)==dg and sum(e[6:])==rd});val=hp(f)*de**dp
 for ee,co in val.lift().dict().items():
  e=tuple(map(int,ee));exp=(e[0],e[1],e[2]+kp,e[3])
  if exp[2]>0:continue
  ans[exp]=ans.get(exp,K(0))+co
ans={e:co for e,co in ans.items() if co};ck('all negative kappa powers cancel',all(e[2]>=0 for e in ans));ck('complete normalized G wall value zero',all(e[3]==1 for e in ans));normal=B(sum((co*W**e[0]*X**e[1] for e,co in ans.items()),B(0)));ck('original-source G kernel unit',normal==-16*k*W);ck('independent source matches array reader',normal==root(ref['complete_G_wall_normal']))
rows=(2,3,6,21,27,25,9);cols=(0,1,2,3,4,5,7);mat=matrix(B,[[root(wall['matrix'][i,j]) for j in cols] for i in rows]);det=mat.det(algorithm='df');ck('exact seven-source determinant',det==B(-QQ(256)/18225*k)*X*W*r**2*q*m**2 or det==B(QQ(256)/18225*k)*X*W*r**2*q*m**2);aug=det*normal;basis=(B(1),W,X,W*X)
def coords(z):
 f=z.lift();return vector(K,[f.monomial_coefficient(t) for t in (A(1),aw,ax,aw*ax)])
NM=matrix(K,4,4,lambda i,j:coords(aug*basis[j])[i]);norm=NM.det();ck('four-root norm exact protected material',norm==16*(QQ(4096)/18225)**4*k**12*r**8*q**4*m**10)
rec={'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'original_full_Q6':True,'base_derivatives_before_restriction':True,'scope':__doc__};save({'record':rec,'checks':checks,'field':K,'complete_G':G,'source_laurent_coefficients':ans,'normal':normal,'four_root_norm':norm},str(p/'independent_subcritical_single_RR_root_wall_source_reader.sobj'),compress=False);(p/'independent_subcritical_single_RR_root_wall_source_reader.json').write_text(json.dumps(rec,indent=2)+'\n');print('PASS',rec,'NORMAL',normal,'NORM',norm.factor(),flush=True)
