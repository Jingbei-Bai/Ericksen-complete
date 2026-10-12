from sage.all import *
from pathlib import Path
from itertools import combinations
import time,resource,os,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'ordinary_single_axis_dimensionless_blocks.sobj'));P=d['ring'];K=P.base_ring();a,r,H,Z=K.gens();u,v=P.gens();A=matrix(P,d['matrix7']);resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
for i in range(A.nrows()):A[i,5]=P(A[i,5]/(v+2))
vb=4*Z*(H+1)/H-2;U=PolynomialRing(K,'u');uu=U.gen();ph=P.hom([uu,U(vb)],U);C=A.apply_map(ph);rows=(7,11,19,25,28,29,30)
def de(B):return B[0,0]*(B[1,1]*B[2,2]-B[1,2]*B[2,1])-B[0,1]*(B[1,0]*B[2,2]-B[1,2]*B[2,0])+B[0,2]*(B[1,0]*B[2,1]-B[1,1]*B[2,0])
items=[]
for rs in combinations(rows,3):
 f=de(C.matrix_from_rows_and_columns(rs,(2,3,5)))
 if f:items.append((rs,f));print('MINOR',rs,'degree_u',f.degree(),flush=True)
print('COUNT',len(items),flush=True)
gg=gcd([f for _,f in items]) if items else U.zero();print('GCD DEG',gg.degree(),flush=True)
KK=FractionField(PolynomialRing(QQ,names=('a','r')));aa,rr=KK.gens();cy=K.hom([aa,rr,KK.one(),-KK.one()],KK);V=PolynomialRing(KK,'u');vv=V.gen()
cyitems=[]
for rs,f in items:
 ff=V({int(e):cy(c) for e,c in f.dict().items()});cyitems.append(ff)
 if len(cyitems)<=10:print('CYCLIC',rs,ff.factor(),flush=True)
print('CYGCD',gcd(cyitems) if cyitems else 0,flush=True)
save({'ring':U,'matrix':C,'vB':vb,'rows':rows,'minors':items,'gcd':gg,'cyclic_minors':cyitems},str(p/'ordinary_single_axis_Bwall.sobj'),compress=False)
(p/'ordinary_single_axis_Bwall.json').write_text(json.dumps({'minors':len(items),'gcd_degree':int(gg.degree()),'seconds':time.monotonic()-st},indent=2)+'\n');print('DONE',time.monotonic()-st,flush=True)
