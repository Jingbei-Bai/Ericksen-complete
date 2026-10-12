"""Whole relative incoming p12/p32 family, full 37 sources and six normals.

p=(1,0,0,0,0,l), RR=(lambda(v)*l,u,v), with the exact base-dependent
lambda from the certified axis kernel. Normals=(p13,p21,p23,p31,eta,tau),
eta=RR1-lambda(v)*l. All l/u/v are retained, including l=0 and u=0.
The full h,z chain rule is evaluated before this family substitution.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic()
g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));hd=load(str(p/'independent_ratio_global_R_derivative_degreefall.sobj'));ed=load(str(p/'independent_ratio_even_derivative_cyclic_direct.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));kc=load(str(p/'independent_ratio_axis_normal_kernel_certificate.sobj'));axis=load(str(p/'independent_ratio_p12_axis_generic_normal.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+900*1024**2,)*2)
K=g['coefficient_field'];a,b,h,z=K.gens();U=PolynomialRing(K,names=('inc12_l','inc12_u','inc12_v'));l,u,v=U.gens()
lam=U({(0,int(e[0]),int(e[1])):K(c) for e,c in kc['lambda'].dict().items()});A=lam*l
O=PolynomialRing(K,names=tuple('in12_p%d'%j for j in range(6))+tuple('in12_R%d'%j for j in range(3)));ov=O.gens()
def cp(f):return O({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
dcache={}
def factor(dd,der=None):
 out={}
 for ee,c in dd.items():
  e=tuple(map(int,ee));ne=e[1:5];sz=sum(ne)
  if sz>1:continue
  if der is not None:
   key=(c,str(der))
   if key not in dcache:dcache[key]=c.derivative(der)
   c=dcache[key]
  if not c:continue
  for j in range(min(e[6],1-sz)+1):
   key=(sum(e),ne+(j,));val=K(c)*binomial(e[6],j)*A**(e[6]-j)*l**e[5]*u**e[7]*v**e[8]
   out[key]=out.get(key,U.zero())+val
 return {e:f for e,f in out.items() if f}
def slot(d,ne,degree):
 tp=degree-d
 if tp<0 or sum(ne)+tp>1:return None
 if tp==1:return 6
 if not sum(ne):return 0
 return 1+ne.index(1)
def one(ff,degree):
 ans=[U.zero() for _ in range(7)]
 for (d,ne),f in ff.items():
  j=slot(d,ne,degree)
  if j is not None:ans[j]+=f
 return vector(U,ans)
def prod2(ff,gg,degree=5):
 ans=[U.zero() for _ in range(7)]
 for (dl,nl),fl in ff.items():
  for (dr,nr),fr in gg.items():
   ne=tuple(nl[j]+nr[j] for j in range(5));j=slot(dl+dr,ne,degree)
   if j is not None:ans[j]+=fl*fr
 return vector(U,ans)
rows=[];records=[];hc={}
def hf(ref):
 if ref not in hc:hc[ref]=factor(hd['polynomial_pool'][ref])
 return hc[ref]
def add(rec,val):
 assert val[0]==0,('family NOT on source',rec,val[0])
 old=axis['matrix'].row(len(rows))
 for j,oj in enumerate((0,1,2,3,5,6)):
  assert val[j+1].subs({l:0})==U({(0,int(e[0]),int(e[1])):K(c) for e,c in old[oj].dict().items()})
 rows.append(list(val[1:]));records.append(rec)
 print(rec,'nonzero',sum(bool(f) for f in val),'sec',round(time.monotonic()-st,1),flush=True)
 save({'coefficient_field':K,'ring':U,'lambda':lam,'normal_positions':('p13','p21','p23','p31','eta','tau'),'source_records':records,'matrix':matrix(U,rows),'scope':__doc__,'qualification':'Relative initial vanishing and first normal matrix only. All tangent walls retained; no no-lift conclusion from generic rank.'},str(p/'independent_ratio_incoming12_generic_normal.sobj'),compress=False)
Q={qi:cp(oq['quartic_source_rows'][qi]) for qi in (0,1,5,6)};src=[cp(f) for f in g['six_full_sources']]
for qi,f in Q.items():add({'kind':'Q','index':qi},one(factor(f.dict()),int(f.total_degree())))
for j,f in enumerate(src):add({'kind':'R','index':j},one(factor(f.dict()),4))
for item in hd['records']:
 val=sum((sgn*prod2(hf(le),hf(ri)) for sgn,le,ri in item['full_signed_product_DAG']),vector(U,[0]*7))
 add({'kind':'H','direction':item['direction'],'index':item['source_position']},val)
hat=[[cp(g['first_gauged_fields'][i,k]) for k in range(9)] for i in range(3)];base=[[cp(g['base_fields'][i,k]) for k in range(2)] for i in range(3)]
for item in ed['items']:
 i=item['direction'];qi=item['Q_source_index'];f=Q[qi];val=vector(U,[0]*7)
 for k in range(9):val+=prod2(factor(hat[i][k].dict()),factor(f.derivative(ov[k]).dict()))
 for k,bp in enumerate((h,z)):val+=prod2(factor(base[i][k].dict()),factor(f.dict(),bp))
 for co,(_,pos,sgn) in zip(item['global_support_coefficient_dicts'],ed['supports'][i]):
  cpoly=O({tuple(map(int,e)):K(c) for e,c in co.items()})*g['source_unit_coefficients'][i]
  val-=sgn*prod2(factor(cpoly.dict()),factor(src[pos].dict()))
 add({'kind':'E','direction':i,'index':qi},val)
(p/'independent_ratio_incoming12_generic_normal.json').write_text(json.dumps({'sources':len(rows),'whole_family_vanishing':True,'axis_specialization_checks':37*6,'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('DONE',time.monotonic()-st,flush=True)
