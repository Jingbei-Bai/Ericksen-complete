"""Independent source replay on all four wall roots and exact unit norm.

Use the already transformed complete Q6 with exact base chain rule,
contrasting the creator's original Q6 coefficient restriction. Evaluate
formal jets as six-element arrays, then compute seven-det by permutation
and the four-dimensional multiplication norm.
"""
from sage.all import *
from pathlib import Path
from itertools import permutations
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();wall=load(str(p/'independent_subcritical_single_RR_conic_wall.sobj'));E=wall['ring'];EW=E.base_ring();K=EW.base_ring();k,r=K.gens();w=E(EW.gen());x=E.gen();q=1+r;m=1+k;g=load(str(p/'independent_partial_D_Newton_gauges_reader.sobj'));ev=load(str(p/'independent_delta_subcritical_even_full_source.sobj'));old=g['ring'].base_ring();R=ev['ring'];OK=R.base_ring();a0,r0,H0,Z0=OK.gens();ph=old.ring().hom([k/r,k*q/(r*m),2/r,-q/r],K);po=OK.ring().hom([k/r,r,2,-1],K);cache={};checks=[]
def ck(label,value):assert bool(value),label;checks.append(label)
def nk(c):
 c=K(c)
 if not c:return K.zero()
 nn,dd=c.numerator(),c.denominator();gg=nn.gcd(dd);nn,dd=nn//gg,dd//gg;lc=dd.leading_coefficient();return K(nn/lc)/K(dd/lc)
def ne(v):
 v=E(v);return E([EW([nk(v[i][j]) for j in range(2)]) for i in range(2)])
def cv(c,field):
 c=field(c);ky=(str(field),c)
 if ky not in cache:
  hh=ph if field==old else po;cache[ky]=E(nk(hh(c.numerator())/hh(c.denominator())))
 return cache[ky]
def add(a,b):return [ne(u+v) for u,v in zip(a,b)]
def mul(a,b):
 out=[E.zero()]*6
 for i in range(6):
  if not a[i]:continue
  ki,ti=divmod(i,2)
  for j in range(6):
   if not b[j]:continue
   kj,tj=divmod(j,2)
   if ki+kj<=2 and ti+tj<=1:out[2*(ki+kj)+ti+tj]=ne(out[2*(ki+kj)+ti+tj]+a[i]*b[j])
 return out
def scale(a,c):return [ne(v*c) for v in a]
one=[E.one(),E.zero(),E.zero(),E.zero(),E.zero(),E.zero()];taujet=[E.zero(),E.one(),E.zero(),E.zero(),E.zero(),E.zero()];kapjet=[E.zero(),E.zero(),E.one(),E.zero(),E.zero(),E.zero()]
ker=wall['kernel'][0];ck('true whole28 wall kernel',wall['matrix']*vector(E,ker)==0);ck('kernel free t and no kappa',ker[6]==1 and ker[7]==0)
ee=[add(scale(one,-2*x*w),scale(taujet,ker[0])),add(scale(one,-w/2),scale(taujet,ker[1])),add(scale(one,2*x),scale(taujet,ker[2]))];alpha=E.one();beta=E(4)/3
images=[scale(one,x),one,scale(one,-x),scale(one,w),add(scale(taujet,-x/beta),scale(mul(kapjet,ee[2]),-1/beta)),add(scale(taujet,x*w/beta),scale(mul(kapjet,add(scale(ee[1],x),scale(ee[0],-1))),1/alpha)),add(scale(taujet,alpha*w/beta),mul(kapjet,ee[1])),taujet,one]
it=g['items'][1];Q=next(it['full_source'] for it in ev['items'] if int(it['Q_index'])==6);v=R.gens();parts=[(it['field_numerators'][j],Q.derivative(v[j])) for j in range(9)]
for j,param in enumerate((H0,Z0)):
 ss=r0 if j==0 else r0/(1+r0);der=R({tuple(map(int,ex)):ss*c.derivative(param) for ex,c in Q.dict().items() if c.derivative(param)});parts.append((it['base_numerators'][j],der))
fc={}
def fmap(f):
 ky=id(f)
 if ky in fc:return fc[ky]
 out={};field=f.parent().base_ring()
 for ex,c in f.dict().items():
  ex=tuple(map(int,ex));dg=sum(ex);rd=sum(ex[6:]);jet=one
  for image,n in zip(images,ex):
   for unused in range(n):jet=mul(jet,image)
  ky=(dg,rd);out[ky]=add(out.get(ky,[E.zero()]*6),scale(jet,cv(c,field)))
 fc[id(f)]=out;return out
ans={}
for left,right in parts:
 for (di,ji),vi in fmap(left).items():
  for (dj,jj),vj in fmap(right).items():
   dg=di+dj;rd=ji+jj;basekap=13-2*dg-rd;dp=6-dg
   if dp<0:raise AssertionError('degree exceeds complete six')
   if dp>1:continue
   jet=mul(vi,vj)
   for power,val in enumerate(jet):
    ki,ti=divmod(power,2);kp=basekap+ki;tp=ti+dp
    if kp>0 or tp>1:continue
    ky=(kp,tp);ans[ky]=ne(ans.get(ky,E.zero())+val*ker[3]**dp)
ans={ky:v for ky,v in ans.items() if v};ck('all negative-kappa full-source coefficients zero',all(ky[0]>=0 for ky in ans));ck('full source value zero on wall',ans.get((0,0),E.zero())==0);Gker=ans.get((0,1),E.zero());ck('independent complete G wall normal',Gker==-16*k*w)
rows=(2,3,6,21,27,25,9);cols=(0,1,2,3,4,5,7);A=wall['matrix'].matrix_from_rows_and_columns(rows,cols);det=E.zero();terms=0
for perm in permutations(range(7)):
 if any(not A[i,perm[i]] for i in range(7)):continue
 term=E.one();sgn=(-1)**sum(perm[i]>perm[j] for i in range(7) for j in range(i+1,7))
 for i in range(7):term=ne(term*A[i,perm[i]])
 det=ne(det+sgn*term);terms+=1
expected=E(QQ(256)/18225)*k*x*w*r*r*q*m*m;ck('true seven determinant protected',det==expected or det==-expected);aug=ne(det*Gker);basis=(E.one(),w,x,w*x)
def coords(v):
 v=E(v);return vector(K,(v[0][0],v[0][1],v[1][0],v[1][1]))
N=matrix(K,4,4,lambda i,j:coords(ne(aug*basis[j]))[i]);norm=nk(N.det());expectednorm=16*(QQ(4096)/18225)**4*k**12*r**8*q**4*m**10;ck('independent four-root multiplication norm',norm==expectednorm)
rel=load(str(p/'independent_subcritical_single_RR_complete_E_relative_dual.sobj'));ck('fullgeneric nonreduced relative inclusion',rel['remainder']=={});print('PASS',len(checks),'seven determinant',det,'Gker',Gker,'norm',norm.factor(),'sec',time.monotonic()-st,flush=True)
save({'ring':E,'checks':checks,'rows':rows,'columns':cols,'seven_det':det,'complete_G_wall_normal':Gker,'augmented_det':aug,'four_root_norm':norm,'source_jet':ans,'scope':__doc__},str(p/'independent_subcritical_single_RR_complete_E_wall_reader.sobj'),compress=False);(p/'independent_subcritical_single_RR_complete_E_wall_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'independent_det_terms':terms,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
