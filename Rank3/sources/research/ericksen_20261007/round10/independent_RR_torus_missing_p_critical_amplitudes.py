"""All original critical H consistency minors on the three missing-p supports."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_RR_torus_missing_p_sources.sobj'));old=d['field'];a0,r0,H0,Z0,x0,b0,c0=old.gens()
K=PolynomialRing(QQ,names=('ca','cr','cH','cZ','cy','cb','cc')).fraction_field();a,r,H,Z,y,b,c=K.gens();BR=K.ring();A=PolynomialRing(QQ,names=('da','dr','dy')).fraction_field();aa,rr,yy=A.gens();P=PolynomialRing(A,names=('amp_b','amp_c'));pb,pc=P.gens();coefph=BR.hom([P(aa),P(rr),P.one(),P.one(),P(yy),pb,pc],P);items=[]
for chart,it in enumerate(d['items']):
 M=it['matrix'];rs=it['generic_pivot_rows'][:3];C=M.matrix_from_rows_and_columns(rs,(0,1,2));det=C.det();adj=C.adjugate();rd=vector(old,[M[i,3] for i in rs]);ro=vector(old,[M[i,4] for i in rs]);qb=vector(old,M[6,:3].list());qt=det*M[6,3]-qb*adj*rd;qz=det*M[6,4]-qb*adj*ro
 xf=[y/b,y,y*b][chart];ph=old.ring().hom([a,r,H,Z,xf,b,c],K)
 def cv(v):v=old(v);return ph(v.numerator())/ph(v.denominator())
 rows=[];sources=[]
 for j in range(7,25):
  hb=vector(old,M[j,:3].list());ht=det*M[j,3]-hb*adj*rd;hz=det*M[j,4]-hb*adj*ro;src=(qt*hz-qz*ht)/det;f=cv(src)
  if not f:continue
  nu=f.numerator();groups={}
  for ex,co in nu.dict().items():
   e=list(map(int,ex));key=tuple(e[2:4]);e[2]=e[3]=0;groups.setdefault(key,{})[tuple(e)]=co
  pol=[(key,coefph(BR(dd))) for key,dd in groups.items()];mins=[min(int(ex[k]) for key,g in pol for ex in g.dict()) for k in range(2)];unit=pb**mins[0]*pc**mins[1];cofs=[]
  for key,g in pol:
   g=g//unit;assert set(tuple(map(int,e)) for e in g.dict())<=set(((0,0),(2,0),(0,2))), (it['label'],j,g)
   row=[g.monomial_coefficient(mon) for mon in (P.one(),pb**2,pc**2)];rows.append(row);cofs.append((key,g))
  sources.append({'H_row':j-7,'original_rows':tuple(rs)+(6,j),'original_determinant':src,'mapped_determinant':f,'removed_amplitude_unit':unit,'base_coefficients':cofs})
 NN=matrix(A,rows);ker=vector(A,[1,-1-rr,-(1+rr)/rr]);assert NN*ker==0
 print('AMPLITUDE',it['label'],NN.dimensions(),'seconds',time.monotonic()-st,flush=True)
 # Normalize protected material units in each row; retain every y factor.
 U=PolynomialRing(QQ,'uy');uy=U.gen();uph=A.ring().hom([U.one(),U.one(),uy],U);conseq=[];records=[]
 for i in range(NN.nrows()):
  for j in range(i+1,NN.nrows()):
   f=NN[i,1]*NN[j,2]-NN[i,2]*NN[j,1]
   if not f:continue
   fac=f.numerator().factor();core=A.ring().one();unitq=A(fac.unit())
   for g,e in fac:
    if g.degree(A.ring().gen(2))>0:core*=g**e
    else:unitq*=A(g)**e
   assert core.degree(A.ring().gen(0))==0 and core.degree(A.ring().gen(1))==0
   pol=uph(core);conseq.append(pol);records.append({'pair':(i,j),'minor':f,'direction_part':pol,'removed_material_unit':unitq/f.denominator()})
   if len(conseq)>=2 and gcd(conseq).degree()<=1:break
  if len(conseq)>=2 and gcd(conseq).degree()<=1:break
 gg=gcd(conseq);print('DIRECTION GCD',it['label'],gg,'pairs',len(conseq),'seconds',time.monotonic()-st,flush=True)
 cofs=U.ideal(conseq).gens() if False else None
 # Repeated xgcd produces explicit original-two-minor cofactors without saturation.
 gcdv=U.zero();mult=[]
 for f in conseq:
  g,ca,cb=gcdv.xgcd(f);mult=[ca*v for v in mult]+[cb];gcdv=g
 assert sum((c*f for c,f in zip(mult,conseq)),U.zero())==gcdv
 items.append({'label':it['label'],'amplitude_matrix':NN,'amplitude_kernel':ker,'sources':sources,'direction_minor_records':records,'direction_gcd':gcdv,'direction_bezout':mult})
 save({'field':A,'amplitude_ring':P,'direction_ring':U,'items':items,'scope':__doc__},str(p/'independent_RR_torus_missing_p_critical_amplitudes.sobj'),compress=False)
print('DONE','seconds',time.monotonic()-st,flush=True)
