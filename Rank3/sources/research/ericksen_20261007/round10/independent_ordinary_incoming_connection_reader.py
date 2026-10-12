"""Independent full original-source reader of ordinary incoming transport.

Only the connection and exact identities are certified. Finiteness of the
ordinary boundary fibers remains a separate hypothesis, not a conclusion.
"""
from sage.all import *
from pathlib import Path
import time,json,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic();checks=[]
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+850*1024**2,)*2)
raw=load(str(p/'independent_ratio_p_boundary_second_gauge.sobj'));ref=load(str(p/'independent_ordinary_incoming_horizontal_fields.sobj'))
F=ref['field'];a,r,H,Z,l,u,v=F.gens();q=1+r;m=1+a*r;K=raw['coefficient_field'];mp=K.hom([a,a*q/m,H/r,q*Z/r],F)
beta=H*(Z-1)/(Z*(H+1));lam=((Z-H)*v-2*H*(Z+1))/(Z*(H+1));S=(Z-H-1)*v-2*H*Z
pv=[F(1),F(0),F(0),F(0),F(0),l,lam*l,u,v]
def ck(label,value):assert value,label;checks.append(label)
def pull(ff,degree):
 out=F(0)
 for ee,c in ff.dict().items():
  e=tuple(map(int,ee))
  if sum(e)!=degree or any(e[j] and not pv[j] for j in range(9)):continue
  out+=mp(c)*prod(pv[j]**e[j] for j in range(9) if e[j])
 return out
items=[]
for i,idx,rpos,sign in [(0,0,3,1),(2,5,1,-1)]:
 ps=raw['first_gauged_fields'][i,0].parent().gens()[idx]
 full=[ps*raw['first_gauged_fields'][i,j] for j in range(9)]
 full[6+i]-=sign*raw['source_unit_coefficients'][i]*raw['six_full_sources'][rpos]
 # Check all top-degree cancellations before restricting to the family.
 for j,f in enumerate(full):ck('global_degree_bound_%d_%d'%(i,j),f.total_degree()<=3)
 ff=[pull(f,3)/pv[idx] for f in full];proj=[ff[j]-pv[j]*ff[0] for j in range(9)]
 bb=[r*pull(ps*raw['base_fields'][i,0],2)/pv[idx],r/q*pull(ps*raw['base_fields'][i,1],2)/pv[idx]]
 old=next(x for x in ref['items'] if x['direction']==i)
 for label,got,want in [('raw',ff,old['raw_fields']),('projective',proj,old['projective_fields']),('base',bb,old['base_fields'])]:
  for j,(x,y) in enumerate(zip(got,want)):ck('source_%d_%s_%d'%(i,label,j),x==y)
 for j in (1,2,3,4,5):ck('p_tangent_constant_l_%d_%d'%(i,j),proj[j]==0)
 ck('RR1_family_tangency_%d'%i,proj[6]==l*(lam.derivative(H)*bb[0]+lam.derivative(Z)*bb[1]+lam.derivative(v)*proj[8]))
 items.append({'direction':i,'fields':proj,'base':bb})
B=matrix(F,[x['base'] for x in items]);A=ref['connection']
ck('true_connection_multiplyback',B*A==matrix(F,[[x['fields'][j] for j in (5,7,8)] for x in items]))
ck('determinant',B.det()==-l*H**2*(v+2)*S/(H+1)**2)
def deriv(f,i):return f.derivative((H,Z)[i])+sum(A[i,j]*f.derivative(x) for j,x in enumerate((l,u,v)))
for j in range(3):ck('flat_%d'%j,deriv(A[1,j],0)==deriv(A[0,j],1))
J=u**2/(beta**2*(v+2)*(v-2*H))
for i in range(2):ck('rational_first_integral_%d'%i,deriv(J,i)==0)
ck('scalar_v_equation',A[1,2]==(v-2*H)*(v+2)/S)
residue=-H/(2*(H+1))
rat=-(H+1)*v/((v-2*H)*(v+2)**2)
ck('rational_differential_partial_fraction',rat==residue/(v-2*H)-residue/(v+2)-1/(v+2)**2)
ck('v_plus_horizontal_H',A[0,2].subs({v:2*H})==2)
ck('v_plus_horizontal_Z',A[1,2].subs({v:2*H})==0)
ck('v_plus_u_exponential',((A[1,1]/u)-beta.derivative(Z)/beta).subs({v:2*H})==-1/(2*H))
singular=[]
for it in items:
 f=S.derivative(H)*it['base'][0]+S.derivative(Z)*it['base'][1]+S.derivative(v)*it['fields'][8]
 wall=f.subs({v:2*H*Z/(Z-H-1)})
 singular.append({'direction':it['direction'],'S_derivative':f,'S_wall_value':wall})
 print('S_WALL_TRANSVERSAL',it['direction'],wall.factor() if wall else 0,flush=True)
out={'field':F,'checks':checks,'items':items,'connection':A,'first_integrals':(l,J),'log_residue':residue,'singular_S':singular,'scope':__doc__}
save(out,str(p/'independent_ordinary_incoming_connection_reader.sobj'),compress=False)
(p/'independent_ordinary_incoming_connection_reader.json').write_text(json.dumps({'status':'PASS','checks':len(checks),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
print('PASS',len(checks),'seconds',time.monotonic()-st,flush=True)
