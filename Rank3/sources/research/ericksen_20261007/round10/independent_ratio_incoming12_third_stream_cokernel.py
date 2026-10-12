"""Stream complete-source third cokernels BEFORE reducing by P(y).

Only the previously saved range pivot Delta is inverted. Coefficient-field
row normalizations are separately recorded, and no pivot is silently a unit.
"""
from sage.all import *
from pathlib import Path
import os,resource,time,json,gc,sys
p=Path(__file__).resolve().parent;st=time.monotonic();c=load(str(p/'independent_ratio_incoming12_generic_nullwall_all_cokernel.sobj'));j=load(str(p/'independent_ratio_incoming12_nullwall_second_compatibility.sobj'));K=c['ring'].base_ring();dd=-c['relation'][0];P=j['primitive_quadratic'];Y=PolynomialRing(K,'stream_y');y=Y.gen();PP=Y(list(P));co={int(it['source_row']):it for it in c['items']};prefix=p/'independent_ratio_incoming12_third_stream_cokernel';src=p/'independent_ratio_incoming12_all_third_raw'
base={rr:load(str(src)+'.row%d.sobj'%rr) for rr in (1,17,5)}
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+900*1024**2,)*2)
zero=(K(0),)*8
def add(x,y):return tuple(x[i]+y[i] for i in range(8))
def smul(u,v):
 out=[]
 for i in range(4):out.extend([u[0]*v[2*i]+dd*u[1]*v[2*i+1],u[0]*v[2*i+1]+u[1]*v[2*i]])
 return tuple(out)
def divide(u,v):
 norm=v[0]**2-dd*v[1]**2
 return ((u[0]*v[0]-dd*u[1]*v[1])/norm,(u[1]*v[0]-u[0]*v[1])/norm)
def atomic(obj,target):
 tmp=Path(str(target)+'.tmp.sobj');save(obj,str(tmp),compress=False);os.replace(tmp,target)
metas=[];selection=list(range(37)) if len(sys.argv)==1 else [int(q) for q in sys.argv[1].split(',')]
for row in selection:
 target=Path(str(prefix)+'.row%d.sobj'%row)
 if target.exists():print('resume',row,flush=True);continue
 it=load(str(src)+'.row%d.sobj'%row);terms=[it,base[1],base[17],base[5]];cf=co[row]['cofactors'];ratios=[divide(v,cf[0]) for v in cf];vals=[]
 for slot in range(6):
  value=zero
  for term,cc in zip(terms,ratios):
   vv=term['unknown_coefficients'][slot] if slot<5 else term['constant'];value=add(value,smul(cc,vv))
  vals.append(value)
 for col in range(6):
  value=zero
  for term,cc in zip(terms,ratios):value=add(value,smul(cc,term['third_range_coefficients'][col]))
  assert value==zero
 # eta2=l*E2 and dl2=l*L2; other correction variables unchanged.
 vals=[smul((0,1),v) if q in (1,2) else v for q,v in enumerate(vals)]
 nz=[v for v in vals if any(v)]
 if not nz:
  pols=[Y.zero()]*6;pivot=K.one();parity=None
 else:
  parity=0 if any(nz[0][2*i] for i in range(4)) else 1
  if not all(not v[2*i+1-parity] for v in vals for i in range(4)):
   print('PARITY',row,[[j for j,cc in enumerate(v) if cc] for v in vals],flush=True)
  assert all(not v[2*i+1-parity] for v in vals for i in range(4))
  pols=[Y([v[2*i+parity] for i in range(4)]) for v in vals]
  pivot=pols[2][0] if pols[2][0] else next(cc for f in pols for cc in f if cc)
  pols=[f/pivot for f in pols]
 divs=[f.quo_rem(PP) for f in pols];reduced=[rr for qq,rr in divs]
 payload={'coefficient_field':K,'ring':Y,'quadratic':PP,'source_row':row,'source_range_ratios':ratios,'parity':parity,'row_pivot':pivot,'unreduced_coefficients':pols,'quotients':[qq for qq,rr in divs],'coefficients':reduced,'unknowns':['k2','Eta2','Dl2','du2','dv2','constant'],'unknown_scaling':'eta2=l*Eta2, dl2=l*Dl2','scope':__doc__}
 atomic(payload,target)
 meta={'row':row,'supports':[q for q,f in enumerate(reduced) if f],'coefficient_termcounts':[sum(len(K(cc).numerator().dict()) for cc in f) for f in reduced],'seconds':time.monotonic()-st};metas.append(meta);print(meta,flush=True)
 del it,terms,payload,pols,reduced,divs,vals;gc.collect()
Path(str(prefix)+'.json').write_text(json.dumps({'processed':metas,'scope':__doc__},indent=2)+'\n');print('DONE',time.monotonic()-st,flush=True)
