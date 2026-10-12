"""Read-only bounded audit of u-content and degree-preserving state slices."""
from pathlib import Path
import time,json
from sage.all import load,save,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic()
print('loading primitive source',flush=True)
d=load(str(out/'pencil_base_content_and_klein.sobj'));R=d['fiber_ring'];K=R.base_ring();u,v,s=R.gens()
print('loaded',round(time.monotonic()-start,2),flush=True)
summary=[];coeffs=[]
for i,p in enumerate(d['primitive_rows']):
 grouped={}
 for mon,c in p.dict().items():
  grouped.setdefault(mon[0],{})[tuple(mon[1:])]=c
 summary.append({'source':i,'u_degree':int(p.degree(u)),'v_degree':int(p.degree(v)),'s_degree':int(p.degree(s)),
  'u_coefficient_support':{str(j):[tuple(m) for m in cc] for j,cc in grouped.items()}})
 coeffs.append(grouped)
print(json.dumps(summary,indent=2),flush=True)
payload={'parameters':d['parameters'],'base_ring':d['base_ring'],'coefficient_dicts_by_u_power':coeffs,
 'source_file':'pencil_base_content_and_klein.sobj','records':summary}
save(payload,str(out/'pencil_state_slice_probe.sobj'))
(out/'pencil_state_slice_probe.json').write_text(json.dumps({'records':summary,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
