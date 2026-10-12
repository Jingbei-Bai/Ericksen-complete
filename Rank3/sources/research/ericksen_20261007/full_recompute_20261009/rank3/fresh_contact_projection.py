"""Fresh fixed-degree projection arithmetic from newly regenerated F/G.

No historical resultants, inverse matrices, Taylor jets, factors, Bezout
coefficients or PASS reports are imported. The finite reduction is a ratio
parameter reduction, justified by the characteristic-zero degree/divisibility
proof below; it is not a material-spectrum test.
"""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;start=time.monotonic();checks=[];results={}
def ck(name,value):
 if not value:raise AssertionError(name)
 checks.append(name)
def progress(name,**kw):print(name,round(time.monotonic()-start,2),kw,flush=True)
data=load(str(out/'fresh_FG_sources.sobj'));src=data['sources'];C=data['ring'];H,Z,S=C.gens()
U=PolynomialRing(QQ,'h');h=U.gen();K=U.fraction_field();SS=PolynomialRing(K,'s');ss=SS.gen()
coeff=[]
for f in src:
 ds={}
 for e,c in f.dict().items():
  he,ze,se=map(int,e);ds[ze,se]=ds.get((ze,se),U.zero())+c*h**he
 coeff.append(ds)
def jets(idx,z0,s0,cap):
 ans={};z0=K(z0);s0=K(s0)
 zp=[z0**j for j in range(max(e[0] for e in coeff[idx])+1)]
 sp=[s0**j for j in range(max(e[1] for e in coeff[idx])+1)]
 for (ze,se),c in coeff[idx].items():
  for i in range(min(ze,cap)+1):
   zpart=c*binomial(ze,i)*zp[ze-i]
   if not zpart:continue
   for j in range(min(se,cap-i)+1):
    val=zpart*binomial(se,j)*sp[se-j]
    if val:ans[i,j]=ans.get((i,j),K.zero())+val
 return {e:c for e,c in ans.items() if c}
order_records=[]
for label,z0,s0,expect,required in [('Z=0',0,0,[4,3,6],12),('Z=H',h,1,[5,5,8],25),('Z=1',1,0,[5,5,8],25),('Z=H+1',h+1,1,[4,5,8],26)]:
 js=[jets(i,z0,s0,8) for i in range(3)]
 orders=[min(map(sum,d),default=9) for d in js]
 ck('finite point orders '+label,orders[:2]==expect[:2] and orders[2]>=expect[2])
 if label!='Z=H+1':ck('intersection bound '+label,min(orders[0]*orders[1],orders[0]*orders[2])>=required)
 else:ck('H+1 contact bound',orders[0]*orders[2]>=required)
 order_records.append(dict(line=label,orders=orders,common_resultant_lower_order=required))
 progress('finite Taylor orders',line=label,orders=orders)

# Resolve the two exceptional tangents at Z=H+1, S=1.
alpha=(2*h+1)/(2*(h*h+h+1));beta=(2*h+1)/(3*h*(h+1));blow=[]
for idx in [0,1]:
 jj=jets(idx,h+1,1,11);mult=[4,5][idx]
 for tangent in [alpha,beta]:
  ans={}
  for (i,j),c in jj.items():
   ep=i+j-mult
   ck('blowup nonnegative '+str((idx,i,j)),ep>=0)
   for nn in range(min(j,6-ep)+1):
    val=c*binomial(j,nn)*tangent**(j-nn)
    if val:ans[ep,nn]=ans.get((ep,nn),K.zero())+val
  ans={e:c for e,c in ans.items() if c};blow.append(ans)
f,g=blow[0],blow[2]
ck('first tangent simple F',f.get((0,1),0)!=0)
slope=-f.get((1,0),0)/f[0,1]
ck('first tangent exact slope',slope==-(h-1)*(4*h*h+7*h+4)/(8*(h*h+h+1)**3))
ck('first tangent cubic cancellation',sum(c*slope**j for (i,j),c in g.items() if i+j==3)==0)
ck('first tangent G order3',min(map(sum,g))==3)
ck('second tangent contribution2',min(map(sum,blow[1]))*min(map(sum,blow[3]))>=2)
ck('two distinct tangents',alpha!=beta)
results['finite_point_orders']=order_records;results['exceptional_slope']=str(slope)
progress('exceptional tangents resolved')

def pval(f,p):
 if not f:return Infinity
 v=0
 while True:
  qq,rr=f.quo_rem(p)
  if rr:return v
  f=qq;v+=1
nonguard=[]
for label,z0,required in [('Z=-1',-1,20),('Z=-H',-h,20),('(H+1)Z=H',h/(h+1),6)]:
 ev=C.hom([SS(h),SS(z0),ss],SS);cache={}
 def zd(fi,n):
  if (fi,n) not in cache:cache[fi,n]=ev(src[fi].derivative(Z,n))
  return cache[fi,n]
 rec=dict(line=label,required=required,intersections=[])
 for oi in [1,2]:
  common=gcd(zd(0,0),zd(oi,0));total=0;factors=[]
  for fac,power in common.factor():
   oo=[]
   for fi in [0,oi]:
    cap=int(pval(zd(fi,0),fac));oo.append(min(nn+pval(zd(fi,nn),fac) for nn in range(cap+1)))
   total+=int(fac.degree())*oo[0]*oo[1];factors.append(dict(factor=str(fac),orders=list(map(int,oo))))
  ck('nonguard order '+label+'/'+str(oi),total>=required)
  rec['intersections'].append(dict(other=oi,lower_order=int(total),factors=factors))
 nonguard.append(rec);progress('nonguard orders',line=label)
results['nonguard']=nonguard

# Reconstruct weighted fixed-size Sylvester initial matrices; their weights
# follow the coefficient-degree pattern, and every inequality is checked.
allm=[];bounds=[]
for oi in [1,2]:
 m=11;n=[0,14,25][oi];size=m+n
 rows=[12-i for i in range(n)]+[(15 if oi==1 else 26)-i for i in range(m)]
 cols=list(range(size));mats=[matrix(K,size,size,0) for unused in range(4)]
 for i in range(size):
  fi=0 if i<n else oi
  for j in range(size):
   se=m-j+i if i<n else n-j+i-n
   for (ze,see),c in coeff[fi].items():
    if see!=se:continue
    ep=rows[i]+cols[j]-ze
    ck('weighted Sylvester support '+str((oi,i,j,ze)),ep>=0)
    if ep<4:mats[ep][i,j]=c
 bounds.append(sum(rows)+sum(cols));allm.append(mats)
 ck('coarse bound '+str(oi),bounds[-1]==[487,861][oi-1])
progress('weighted Sylvester matrices rebuilt')
A0=allm[0][0];cols=list(A0.pivots());pivot_rows=list(A0.matrix_from_columns(cols).transpose().pivots())
free=[i for i in range(25) if i not in cols];rem=[i for i in range(25) if i not in pivot_rows]
ck('infinity initial rank23',len(cols)==23 and len(pivot_rows)==23)
B=A0.matrix_from_rows_and_columns(pivot_rows,cols);BI=B.inverse();lifts=[];schurs=[]
for n in range(4):
 xx=matrix(K,25,2,0)
 if n==0:
  for j,fcol in enumerate(free):xx[fcol,j]=1
 rhs=sum((allm[0][j]*lifts[n-j] for j in range(1,n+1)),matrix(K,25,2,0))
 correction=-BI*(rhs+A0*xx).matrix_from_rows(pivot_rows)
 for i,idx in enumerate(cols):xx[idx,:]=correction[i,:]
 residual=A0*xx+rhs
 ck('Schur pivot vanishes '+str(n),residual.matrix_from_rows(pivot_rows)==0)
 lifts.append(xx);schurs.append(residual.matrix_from_rows(rem))
ck('Schur order0',schurs[0]==0);ck('Schur order1',schurs[1]==0)
A2,A3=schurs[2:4]
ck('Schur determinant order4 zero',A2.det()==0)
ck('Schur determinant order5 zero',A2[0,0]*A3[1,1]+A3[0,0]*A2[1,1]-A2[0,1]*A3[1,0]-A3[0,1]*A2[1,0]==0)
ck('strict degree bound481',bounds[0]-6==481)
save(dict(pivot_rows=pivot_rows,pivot_columns=cols,constant_inverse=BI,kernel_lifts=lifts,small_schur_coefficients=schurs,series=allm),str(out/'fresh_infinity_schur.sobj'),compress=False)
results['infinity']=dict(initial_rank=23,pivot_rows=pivot_rows,pivot_columns=cols,coarse_bound=487,vanishing_order=6,strict_bound=481)
progress('fresh infinity Schur bound481')

# Exact modular resultant computation by polynomial interpolation of explicitly
# padded determinants. This does not load or compare to saved resultants.
prime=32003;hv=3;Ff=GF(prime);UU=PolynomialRing(Ff,'zz');zz=UU.gen();modco=[]
for fi in range(3):
 ds=[{} for j in range([11,14,25][fi]+1)]
 for (he,ze,se),c in src[fi].dict().items():
  ck('good coefficient denominator',int(c.denominator())%prime!=0)
  ds[int(se)][int(ze)]=ds[int(se)].get(int(ze),Ff.zero())+Ff(c)*Ff(hv)**int(he)
 modco.append([UU(q) for q in ds]);ck('S degree preserved '+str(fi),modco[-1][-1]!=0)
rs=[]
for oi,bound in [(1,487),(2,861)]:
 m=11;n=[0,14,25][oi];size=m+n;values=[]
 for zvalue in range(bound+1):
  aa=[f(zvalue) for f in modco[0]];bb=[f(zvalue) for f in modco[oi]];mat=matrix(Ff,size,size,0)
  for i in range(n):
   for j in range(m+1):mat[i,i+j]=aa[m-j]
  for i in range(m):
   for j in range(n+1):mat[n+i,i+j]=bb[n-j]
  values.append((Ff(zvalue),mat.det()))
 rp=UU.lagrange_polynomial(values);rs.append(rp)
 progress('fresh modular padded resultant',other=oi,degree=int(rp.degree()))
D0=zz**12*(zz+1)**20*(zz+hv)**20*(zz-hv)**25*(zz-1)**25*(zz-hv-1)**26*((hv+1)*zz-hv)**6
quot=[]
for i,rp in enumerate(rs):
 q,rem0=rp.quo_rem(D0);ck('exact modular D0 division '+str(i),rem0==0);quot.append(q)
g,bez1,bez2=xgcd(*quot);ck('residual Bezout unit',g==1);ck('fresh Bezout multiplication',bez1*quot[0]+bez2*quot[1]==1)
ck('first degree preserved',rs[0].degree()==481 and quot[0].degree()==347)
save(dict(prime=prime,H=hv,resultants=rs,D0=D0,quotients=quot,bezout=[bez1,bez2]),str(out/'fresh_modular_projection.sobj'),compress=False)
results['modular']=dict(prime=prime,H=hv,determinant_evaluations=[488,862],resultant_degrees=[int(f.degree()) for f in rs],quotient_degrees=[int(f.degree()) for f in quot],gcd=1)
report=dict(status='PASS',checks=checks,check_count=len(checks),results=results,seconds=time.monotonic()-start,no_material_specialization=True,
 scope='Fresh finite orders, exceptional tangents, infinity Schur bound, exact finite-field resultants and Bezout; F/G input was regenerated from explicit manuscript matrices. Source-to-matrix provenance remains separate.')
(out/'fresh_contact_projection.json').write_text(json.dumps(report,indent=2));progress('PASS',count=len(checks))
