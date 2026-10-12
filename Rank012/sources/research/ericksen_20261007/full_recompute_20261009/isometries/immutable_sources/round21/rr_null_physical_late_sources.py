"""Protected Cartesian third-jet elimination BEFORE first-jet substitution.

All polynomial arithmetic is over QQ(i); no intermediate fraction-field matrix
expansion. The final kernel is normalized on its last three columns. Its only
possible intrinsic denominator is det(the first twelve rows), nonzero on the
whole main-direction chart by the independently proved cubic-symbol lemma.
"""
from sage.all import *
from pathlib import Path
import signal,time,json,resource,sys
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(150)
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1400*1024**2,vm+1400*1024**2))
delta=int(sys.argv[1]) if len(sys.argv)>1 else 1;tag='plus' if delta==1 else 'minus'
d=load(str(O/f'rr_null_physical_protected_graph_{tag}.sobj'))
K=d['field'];PR=K.ring();QI=PR.base_ring();ii=QI.gen();T=d['ring']
names=PR.variable_names()+('alpha','beta')+tuple('h'+str(j) for j in range(9))+tuple('z'+str(j) for j in range(12))
P=PolynomialRing(QI,names=names);vv=P.gens();s,t,p,q,c=vv[:5];aa,bb=vv[5:7];hh=vv[7:16];zz=vv[16:]
F=matrix(P,[[s,0,0],[0,t,0],[p,q,c]])
adjF=matrix(P,[[t*c,0,0],[0,s*c,0],[-t*p,-s*q,s*t]])
J=s*t*c;assert adjF*F==J*identity_matrix(P,3)
C=F.transpose()*F;I1=C.trace()
As=vector(P,[aa,ii*aa,0]);At=vector(P,[bb,delta*ii*bb,0])
ee=[vector(P,[int(i==j) for i in range(3)]) for j in range(3)]
H={(0,0):vector(P,hh[:3]),(0,1):vector(P,hh[3:6]),(1,1):vector(P,hh[6:])};H[1,0]=H[0,1]
for a in range(3):H[a,2]=H[2,a]=At.cross_product(F.column(a))-F*As.cross_product(ee[a])
third={}
for ix,ijk in enumerate([(0,0,0),(0,0,1),(0,1,1),(1,1,1)]):third[ijk]=vector(P,zz[3*ix:3*ix+3])
for a in range(3):
 for b in range(a,3):
  va=As.cross_product(ee[a]);vb=As.cross_product(ee[b])
  third[tuple(sorted((a,b,2)))]=At.cross_product(H[a,b])-sum((va[j]*H[j,b]+vb[j]*H[a,j] for j in range(3)),vector(P,[0,0,0]))
def H3(a,b,c):return third[tuple(sorted((a,b,c)))]
DF=[matrix(P,3,3,lambda i,j:H[j,a][i]) for a in range(3)]
DC=[f.transpose()*F+F.transpose()*f for f in DF]
DDF={};DDC={};DDB={}
for a in range(3):
 for b in range(a,3):
  DDF[a,b]=matrix(P,3,3,lambda i,j:H3(j,a,b)[i])
  dd=DDF[a,b].transpose()*F+DF[a].transpose()*DF[b]+DF[b].transpose()*DF[a]+F.transpose()*DDF[a,b]
  DDC[a,b]=dd
  DDB[a,b]=dd*C+C*dd+DC[a]*DC[b]+DC[b]*DC[a]-I1*dd
def mat(dd,a,b):return dd[tuple(sorted((a,b)))]
def spectral(a,b):
 cc=mat(DDC,a,b);ff=mat(DDF,a,b)
 return [cc.trace(),(C*cc).trace()+(DC[a]*DC[b]).trace(),J*(adjF*ff).trace()-(adjF*DF[a]*adjF*DF[b]).trace()]
eq=[];enames=[]
for a,b in [(0,0),(0,1),(1,1)]:eq+=spectral(a,b);enames += [f'I1_{a}{b}',f'I2_{a}{b}',f'J2_{a}{b}']
for label,dd in [('C',DDC),('J2B',DDB)]:
 for i,j in [(0,1),(0,2),(1,2)]:
  eq.append(sum(mat(dd,i,a)[a,j]-mat(dd,j,a)[a,i] for a in range(3)))
  enames.append(f'{label}curl_{i}{j}')
R=PolynomialRing(K,names=('alpha','beta')+tuple('h'+str(j) for j in range(9)))
topr=P.hom(list(PR.gens())+[PR(0)]*23,PR)
tor=P.hom(list(map(R,K.gens()))+list(R.gens())+[R(0)]*12,R)
M=matrix(PR,[[topr(f.derivative(z)) for z in zz] for f in eq]);off=vector(R,[tor(f) for f in eq])
checks={'third_symbol_no_H_dependency':all(f.derivative(z)==P(t) for f,row in zip(eq,M.rows()) for z,t in zip(zz,row))}
assert all(checks.values()),checks
save({'field':K,'ring':R,'flat_ring':P,'equations':eq,'names':enames,'matrix':M,'offset':off,'checks':checks,'delta':delta},str(O/f'rr_null_physical_late_full_source_{tag}.sobj'))
print('late source saved',time.time()-st,flush=True)
kp=O/'rr_null_physical_protected_third_symbol.sobj'
if kp.exists():
 kd=load(str(kp));ker=kd['kernel'];assert kd['matrix']==M
else:
 rows=[list(M.row(i))+[PR(i==j) for j in range(15)] for i in range(15)]
 active=list(range(12));cols=list(range(12));pivs=[]
 for h in range(12):
  choices=[((rows[i][j].total_degree(),len(rows[i][j].dict())),i,j) for i in active for j in cols if rows[i][j]]
  _,i,j=min(choices);active.remove(i);cols.remove(j);p0=rows[i][j];pivs.append((i,j,p0))
  print('pivot',h,'deg',p0.total_degree(),'terms',len(p0.dict()),'sec',time.time()-st,flush=True)
  for l in active+[12,13,14]:
   q0=rows[l][j]
   if not q0:continue
   gg=p0.gcd(q0);a0=p0//gg;b0=q0//gg
   nr=[a0*rows[l][m0]-b0*rows[i][m0] for m0 in range(27)]
   gg=gcd([z for z in nr if z]);rows[l]=[z//gg for z in nr] if gg else nr
  assert all(rows[l][j]==0 for l in active+[12,13,14])
 def cf(x):
  x=K(x);n=x.numerator();de=x.denominator();lc=de.leading_coefficient();return K(n/lc)/K(de/lc)
 ker=matrix(K,[[cf(K(z)/K(rows[i][12+i])) for z in rows[i][12:]] for i in [12,13,14]])
 assert ker*M==0 and ker[:,12:]==identity_matrix(K,3)
 save({'matrix':M,'kernel':ker,'pivots':pivs,'checks':{'kernel':True,'last_three_identity':True}},str(kp))
checks['protected_normalized_kernel']=ker*M==0 and ker[:,12:]==identity_matrix(K,3)
raw=list(ker*off)
save({'field':K,'ring':R,'quadratics_raw_H':raw,'kernel':ker,'checks':checks,'delta':delta},str(O/f'rr_null_physical_before_graph_quadratics_{tag}.sobj'))
print('raw H quadratics saved',time.time()-st,flush=True)
ph=R.hom([T.gen(0),T.gen(1)]+d['graph'],T)
Qraw=[ph(f) for f in raw]
out={'field':K,'ring':T,'raw_quadratics':Qraw,'quadratics':Qraw,'multipliers':[K(1)]*3,'kernel':ker,'checks':checks,'delta':delta}
save(out,str(O/f'rr_null_physical_quadratics_{tag}.sobj'))
rep={'status':'EXACT_WHOLE_MAIN_NULL_INFINITY_QUADRATICS_UNSCALED','checks':checks,'delta':delta,
 'coefficient_terms':[sum(len(z.numerator().dict()) for z in f.coefficients()) for f in Qraw],
 'denominators':sorted(set(str(z.denominator().factor()) for f in Qraw for z in f.coefficients())), 'seconds':time.time()-st}
(O/f'rr_null_physical_quadratics_{tag}.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
