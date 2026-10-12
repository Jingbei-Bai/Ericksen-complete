"""Six globally genuine Q0 derivative degree-fall sources and family normals.

If Rj_top=sign*p_s*Z_sigma(i), set
G_ij=p_s*hatFi(Q0)-sign*u_i*(partial_RRi Q0)*Rj.
The ordinary degree-six layer cancels identically on the full state ring.
No state inverse is introduced, including when the selected p_s is zero.
The six degree<=5 sources are evaluated on the entire incoming12 family.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json
p=Path(__file__).resolve().parent;st=time.monotonic();g=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));oq=load(str(p.parent/'round8/rank3_joint_fiber_leading.sobj'));inc=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'))
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+600*1024**2,)*2)
K=g['coefficient_field'];a,b,h,z=K.gens();U=inc['ring'];l,u,v=U.gens();A=inc['lambda']*l
O=PolynomialRing(K,names=tuple('q0act_p%d'%j for j in range(6))+tuple('q0act_R%d'%j for j in range(3)));ov=O.gens()
def cp(f):return O({tuple(map(int,e)):K(c) for e,c in f.dict().items()})
Q=cp(oq['quartic_source_rows'][0]);src=[cp(f) for f in g['six_full_sources']];hat=[[cp(g['first_gauged_fields'][i,k]) for k in range(9)] for i in range(3)];base=[[cp(g['base_fields'][i,k]) for k in range(2)] for i in range(3)]
def layer(f,d):return O({tuple(e):c for e,c in f.dict().items() if sum(e)==d})
def factor(f):
 out={}
 for ee,c in f.dict().items():
  e=tuple(map(int,ee));ne=e[1:5];sz=sum(ne)
  if sz>1:continue
  for j in range(min(e[6],1-sz)+1):
   key=(sum(e),ne+(j,));out[key]=out.get(key,U.zero())+K(c)*binomial(e[6],j)*A**(e[6]-j)*l**e[5]*u**e[7]*v**e[8]
 return {e:f for e,f in out.items() if f}
def prod2(ff,gg):
 out=[U.zero() for _ in range(7)]
 for (dl,nl),fl in ff.items():
  for (dr,nr),fr in gg.items():
   ne=tuple(nl[j]+nr[j] for j in range(5));tp=5-dl-dr
   if tp<0 or sum(ne)+tp>1:continue
   j=6 if tp else (1+ne.index(1) if sum(ne) else 0)
   out[j]+=fl*fr
 return vector(U,out)
pool={};cache={}
def ref(name,f):
 if name not in pool:pool[name]={tuple(e):c for e,c in f.dict().items()};cache[name]=factor(f)
 return name
qd=[Q.derivative(ov[j]) for j in range(9)];qb=[O({tuple(e):c.derivative(bp) for e,c in Q.dict().items() if c.derivative(bp)}) for bp in (h,z)]
maps=[(0,0,3,1),(0,1,4,1),(1,2,2,1),(1,3,5,-1),(2,4,0,1),(2,5,1,-1)]
items=[]
for i,pp,rpos,sgn in maps:
 ui=K(g['source_unit_coefficients'][i]);ps=ov[pp]
 assert ps*layer(hat[i][6+i],3)==sgn*ui*layer(src[rpos],4)
 assert all(layer(hat[i][j],3)==0 for j in range(9) if j!=6+i)
 dag=[]
 for j in range(9):dag.append((1,ref('p%dF%d_%d'%(pp,i,j),ps*hat[i][j]),ref('Qd%d'%j,qd[j])))
 for j in range(2):dag.append((1,ref('p%dB%d_%d'%(pp,i,j),ps*base[i][j]),ref('Qb%d'%j,qb[j])))
 dag.append((-sgn,ref('u%dQd%d'%(i,i),ui*qd[6+i]),ref('R%d'%rpos,src[rpos])))
 value=sum((sign*prod2(cache[aa],cache[bb]) for sign,aa,bb in dag),vector(U,[0]*7))
 items.append({'direction':i,'p_column':pp,'R_position':rpos,'sign':sgn,'full_signed_product_DAG':dag,'family_value_and_normal':value,'global_degree6_cancellation':True})
 print(i,pp,rpos,'family',len(value[0].dict()),'normal counts',[len(f.dict()) for f in value[1:]],'sec',round(time.monotonic()-st,1),flush=True)
 save({'coefficient_field':K,'state_ring':O,'family_ring':U,'lambda':inc['lambda'],'polynomial_pool':pool,'items':items,'scope':__doc__},str(p/'independent_ratio_global_Q0_active_degreefall.sobj'),compress=False)
(p/'independent_ratio_global_Q0_active_degreefall.json').write_text(json.dumps({'sources':len(items),'family_initial_nonzero':[j for j,it in enumerate(items) if it['family_value_and_normal'][0]],'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('DONE',time.monotonic()-st,flush=True)
