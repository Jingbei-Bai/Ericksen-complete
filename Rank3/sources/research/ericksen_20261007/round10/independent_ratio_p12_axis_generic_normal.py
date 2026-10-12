"""Full relative seven-normal first jets at the p12 axis, all base parameters.

p12=1; tangents are RR2=u, RR3=v. Normals are the five other p entries,
RR1 and tau. All37 full original sources are used, with no state inverse.
No K(u,v) matrix solve, determinant, or rank-specialization inference.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic()
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));hd=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));ed=load(str(p/'independent_ratio_even_derivative_cyclic_direct.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+768*1024**2,)*2)
K=g['coefficient_field'];a,b,h,z=K.gens();U=PolynomialRing(K,names=('relax_u','relax_v'));u,v=U.gens();norm=(1,2,3,4,5,6)
O=PolynomialRing(K,names=tuple('relax_p%d'%j for j in range(6))+tuple('relax_R%d'%j for j in range(3)));ov=O.gens()
def cp(f):return O({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
dcache={}
def factor(dd,der=None):
 out={}
 for ee,c in dd.items():
  e=tuple(map(int,ee));ne=tuple(e[j] for j in norm)
  if sum(ne)>1:continue
  if der is not None:
   key=(c,str(der))
   if key not in dcache:dcache[key]=c.derivative(der)
   c=dcache[key]
  if c:
   key=(sum(e),ne);out[key]=out.get(key,U.zero())+K(c)*u**e[7]*v**e[8]
 return {e:f for e,f in out.items() if f}
def monomial_slot(d,ne,degree):
 tp=degree-d
 if tp<0 or sum(ne)+tp>1:return None
 if tp==1:return 7
 if not sum(ne):return 0
 return 1+ne.index(1)
def one(ff,degree):
 ans=[U.zero() for _ in range(8)]
 for (d,ne),f in ff.items():
  j=monomial_slot(d,ne,degree)
  if j is not None:ans[j]+=f
 return vector(U,ans)
def prod(ff,gg,degree=5):
 ans=[U.zero() for _ in range(8)]
 for (dl,nl),fl in ff.items():
  for (dr,nr),fr in gg.items():
   ne=tuple(nl[j]+nr[j] for j in range(6));slot=monomial_slot(dl+dr,ne,degree)
   if slot is not None:ans[slot]+=fl*fr
 return vector(U,ans)
rows=[];records=[];hc={}
def hf(ref):
 if ref not in hc:hc[ref]=factor(hd['polynomial_pool'][ref])
 return hc[ref]
def add(rec,value):
 assert value[0]==0,'relative boundary value unexpectedly nonzero'
 rows.append(list(value[1:]));records.append(rec)
 print(rec,'nonzero',sum(bool(q) for q in value),'time',round(time.monotonic()-st,1),flush=True)
 save({'coefficient_field':K,'ring':U,'normal_positions':norm+('tau',),'source_records':records,'matrix':matrix(U,rows),'scope':__doc__,'qualification':'All tangent u/v and material/base factors are retained. No ordinary boundary no-lift conclusion follows from this payload alone.'},str(p/'independent_ratio_p12_axis_generic_normal.sobj'),compress=False)
Q={qi:cp(oq['quartic_source_rows'][qi]) for qi in (0,1,5,6)};src=[cp(f) for f in g['six_full_sources']]
for qi,f in Q.items():add({'kind':'Q','index':qi},one(factor(f.dict()),int(f.total_degree())))
for j,f in enumerate(src):add({'kind':'R','index':j},one(factor(f.dict()),4))
for item in hd['records']:
 val=sum((sgn*prod(hf(le),hf(ri)) for sgn,le,ri in item['full_signed_product_DAG']),vector(U,[0]*8))
 add({'kind':'H','direction':item['direction'],'index':item['source_position']},val)
hat=[[cp(g['first_gauged_fields'][i,k]) for k in range(9)] for i in range(3)];base=[[cp(g['base_fields'][i,k]) for k in range(2)] for i in range(3)]
for item in ed['items']:
 i=item['direction'];qi=item['Q_source_index'];f=Q[qi];val=vector(U,[0]*8)
 for k in range(9):val+=prod(factor(hat[i][k].dict()),factor(f.derivative(ov[k]).dict()))
 for k,bp in enumerate((h,z)):val+=prod(factor(base[i][k].dict()),factor(f.dict(),bp))
 for co,(_,pos,sgn) in zip(item['global_support_coefficient_dicts'],ed['supports'][i]):
  cpoly=O({tuple(map(int,e)):K(c) for e,c in co.items()})*g['source_unit_coefficients'][i]
  val-=sgn*prod(factor(cpoly.dict()),factor(src[pos].dict()))
 add({'kind':'E','direction':i,'index':qi},val)
(p/'independent_ratio_p12_axis_generic_normal.json').write_text(json.dumps({'sources':len(rows),'nonzero_counts':[sum(bool(f) for f in row) for row in rows],'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('DONE',time.monotonic()-st,flush=True)
