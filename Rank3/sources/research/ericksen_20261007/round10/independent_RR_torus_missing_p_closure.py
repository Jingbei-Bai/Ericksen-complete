"""Full critical/supercritical source consequences on all three four-p supports."""
from sage.all import *
from pathlib import Path
import time
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_RR_torus_missing_p_super_normal.sobj'));sd=load(str(p/'independent_RR_torus_missing_p_sources.sobj'));K=d['field'];a,r,H,Z,y=K.gens();BR=K.ring();U=PolynomialRing(QQ,'uy');uy=U.gen();items=[];checks=[]
def ck(n,v):assert v,n;checks.append(n)
def consequences(f):
 ck('nonzero true source determinant',f!=0);core=BR.one();un=K(f.numerator().factor().unit());allowed=[r,r+1,H,Z,H+1,Z-H,Z-1,H*Z-H+Z,y]
 def const(v):return v.numerator().total_degree()==0 and v.denominator().total_degree()==0
 for g,e in f.numerator().factor():
  if any(const(K(g)/v) for v in allowed):un*=K(g)**e
  else:core*=g**e
 ck('complete determinant factor replay',K(core)*un/f.denominator()==f)
 for g,e in f.denominator().factor():ck('protected full source denominator',any(const(K(g)/v) for v in allowed))
 groups={}
 for ex,co in core.dict().items():
  ee=tuple(map(int,ex));ck('all spectra retained in residual',ee[0]==ee[1]==0);key=ee[2:4];groups[key]=groups.get(key,U.zero())+QQ(co)*uy**ee[4]
 return {'determinant':f,'core':core,'removed_unit':un/f.denominator(),'base_coefficients':groups}
def bezout(gens):
 cur=U.zero();co=[]
 for f in gens:
  g,ca,cb=cur.xgcd(f);co=[ca*v for v in co]+[cb];cur=g
 ck('original source coefficient Bezout',sum((c*f for c,f in zip(co,gens)),U.zero())==cur)
 return cur,co
for ix,it in enumerate(d['items']):
 N=it['matrix'];rs=sd['items'][ix]['generic_pivot_rows'][:3];C=N.matrix_from_rows_and_columns(rs,(2,3,4));det=C.det();co=vector(K,[N[i,5] for i in rs]);qb=vector(K,N[6,2:5].list());crit=det*N[6,5]-qb*C.adjugate()*co
 cr=consequences(crit);cg,cc=bezout(list(cr['base_coefficients'].values()));print('CRITICAL',it['label'],'factor',crit.factor(),'directiongcd',cg,'seconds',time.monotonic()-st,flush=True)
 pr=it['pivot_rows'];ir=pr[:2];T=N.matrix_from_rows_and_columns(ir,(0,1));td=T.det();V=N[:,2:]-N[:,:2]*T.adjugate()*N.matrix_from_rows_and_columns(ir,(2,3,4,5))/td
 br=pr[2:5];C=V.matrix_from_rows_and_columns(br,(0,1,2));det=C.det();adj=C.adjugate();offset=vector(K,[V[i,3] for i in br]);sr=[];sg=[]
 for j in range(N.nrows()):
  if j in pr[:5]:continue
  f=td*(det*V[j,3]-vector(K,V[j,:3].list())*adj*offset)
  if not f:continue
  rec=consequences(f);rec['source_rows']=tuple(pr[:5])+(j,);sr.append(rec);sg.extend(rec['base_coefficients'].values())
 gg,bc=bezout(sg);print('SUPER',it['label'],'gcd',gg,'source minors',len(sr),'seconds',time.monotonic()-st,flush=True)
 ck('whole critical direction cover',cg==1);ck('whole super direction cover',gg==1)
 items.append({'label':it['label'],'critical_source_rows':tuple(rs)+(6,),'critical':cr,'critical_gcd':cg,'critical_bezout':cc,'super_transverse_rows':ir,'super_transverse_determinant':td,'super_reduced_matrix':V,'super_minors':sr,'super_coefficients':sg,'super_gcd':gg,'super_bezout':bc})
 save({'field':K,'direction_ring':U,'items':items,'checks':checks,'scope':__doc__},str(p/'independent_RR_torus_missing_p_closure.sobj'),compress=False)
print('DONE',len(checks),'seconds',time.monotonic()-st,flush=True)
