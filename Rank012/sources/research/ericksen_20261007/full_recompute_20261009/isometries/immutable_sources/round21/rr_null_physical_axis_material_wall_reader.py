"""Resolve the finite degree-23 real material guard wall exactly."""
from sage.all import *
from pathlib import Path
import time,signal,json,resource,sys
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(100)
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
tag=sys.argv[1] if len(sys.argv)>1 else 'plus';d=load(str(O/f'rr_null_physical_axis_symmetric_resultant_{tag}.sobj'))
Q=d['ring'];z,y=Q.gens();Y=d['resultant'].parent();yp=Y.gen();f=max(d['factorization'],key=lambda x:x[0].degree())[0]
K=NumberField(f,'v');v=K.gen();Z=PolynomialRing(K,'z');zz=Z.gen()
def cv(poly):return sum(K(co)*zz**i*v**j for (i,j),co in poly.dict().items())
fr=cv(d['real']);fs=cv(d['imag_core']);gg=fr.gcd(fs).monic()
checks={'defining_polynomial_irreducible':f.is_irreducible(),'linear_common_factor':gg.degree()==1,'real_source_divisible':fr%gg==0,'imag_source_divisible':fs%gg==0}
assert all(checks.values()),checks
zv=-gg[0]
roots=[r for r,m in f.roots(AA) if r>1 or r< -1]
rows=[];valid=[]
def box(x):
 n=floor(x*10**6);return (QQ(n)/10**6,QQ(n+1)/10**6)
for r in roots:
 hom=K.hom([r],AA);zr=hom(zv);lower=AA(0) if r>1 else -4*r;upper=(r-1)**2
 inside=bool(zr>lower and zr<upper)
 ylo,yhi=box(r);zlo,zhi=box(zr)
 rational_exclusion=bool(zhi<0 or (r>1 and zlo>(yhi-1)**2) or (r< -1 and zhi< -4*yhi))
 checks[f'root_{len(rows)}_rational_interval_exclusion']=rational_exclusion
 rows.append({'y':str(r),'z':str(zr),'lower':str(lower),'upper':str(upper),'inside':inside,'sheet':'same' if r>1 else 'opposite',
              'y_interval':[str(ylo),str(yhi)],'z_interval':[str(zlo),str(zhi)],'rational_interval_exclusion':rational_exclusion})
 if inside:valid.append((r,zr))
checks['all_relevant_real_conjugates_excluded']=not valid
assert all(checks.values()),checks
rep={'status':'EXACT_FINITE_MATERIAL_WALL_AUDIT','checks':checks,'candidate_count':len(roots),'physical_count':len(valid),'candidates':rows,'seconds':time.time()-st}
save({'number_field':K,'gcd':gg,'z_value':zv,'candidates':[(r,K.hom([r],AA)(zv)) for r in roots],'physical_candidates':valid,'checks':checks},str(O/f'rr_null_physical_axis_material_wall_{tag}.sobj'))
(O/f'rr_null_physical_axis_material_wall_{tag}.json').write_text(json.dumps(rep,indent=2),encoding='utf8');print(json.dumps(rep,indent=2),flush=True)
