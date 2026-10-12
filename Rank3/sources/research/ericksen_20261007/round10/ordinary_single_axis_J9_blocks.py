"""Small exact source blocks of J_axis9, u=0, retaining all v."""
from sage.all import *
from pathlib import Path
from itertools import combinations
import resource,os,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'ordinary_single_axis_J9.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+500*1024**2,)*2)
U=d['ring'];K=U.base_ring();oa,ob,oh,oz=K.gens();R=PolynomialRing(K,'v');v=R.gen()
def cv(f):return R({int(e[1]):c for e,c in U(f).dict().items() if not e[0]})
A=d['matrix_u_zero'].apply_map(cv);checks=[]
def ck(n,z):assert bool(z),n;checks.append(n)
def det3(A):return A[0,0]*(A[1,1]*A[2,2]-A[1,2]*A[2,1])-A[0,1]*(A[1,0]*A[2,2]-A[1,2]*A[2,0])+A[0,2]*(A[1,0]*A[2,1]-A[1,1]*A[2,0])
blocks=[]
for label,rs,cs in [('quadratic',(1,2,3,34,35,36),(1,6,8)),('tau',(7,11,25,37),(3,5,7))]:
 items=[]
 for rr in combinations(rs,3):
  f=det3(A.matrix_from_rows_and_columns(rr,cs))
  if not f:continue
  print(label,rr,'degree',f.degree(),'coeffs',len(f.dict()),flush=True)
  items.append({'rows':rr,'polynomial':f})
  # save immediately before adaptive gcd/factor.
  save({'ring':R,'matrix':A,'blocks':blocks,'current':{'name':label,'rows':rs,'columns':cs,'minors':items}},str(p/'ordinary_single_axis_J9_blocks.partial.sobj'),compress=False)
 gg=gcd([it['polynomial'] for it in items]) if items else R.zero()
 block={'name':label,'rows':rs,'columns':cs,'minors':items,'gcd':gg}
 blocks.append(block);print('GCD',label,str(gg),flush=True)
 for e,c in gg.dict().items():print('GCD FACTOR',e,c.factor(),flush=True)
B=PolynomialRing(QQ,names=('a','r'));KK=B.fraction_field();a,r=KK.gens();q=1+r;m=1+a*r;ph=K.hom([a,a*q/m,1/r,-q/r],KK);S=PolynomialRing(KK,'v');vv=S.gen()
def cyc(f):return S({int(e):ph(c) for e,c in R(f).dict().items()})
for block in blocks:
 cf=[cyc(it['polynomial']) for it in block['minors']];block['cyclic_minors']=cf
 cg=gcd(cf);block['cyclic_gcd']=cg
 print('CYCLIC GCD',block['name'],cg,flush=True)
 for it,f in zip(block['minors'],cf):print('CYCLIC FACTOR',block['name'],it['rows'],f.factor(),flush=True)
for ij in [(13,0),(19,2),(0,4)]:print('SINGLE',ij,cyc(A[ij]).factor(),flush=True)
save({'ring':R,'matrix':A,'blocks':blocks,'singleton_entries':{ij:A[ij] for ij in [(13,0),(19,2),(0,4)]},'cyclic_ring':S,'checks':checks},str(p/'ordinary_single_axis_J9_blocks.sobj'),compress=False)
rep={'status':'EXACT BLOCKS','seconds':time.monotonic()-st,'blocks':[{'name':b['name'],'gcd':str(b['gcd']),'cyclic_gcd':str(b['cyclic_gcd'])} for b in blocks]};(p/'ordinary_single_axis_J9_blocks.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep,flush=True)
