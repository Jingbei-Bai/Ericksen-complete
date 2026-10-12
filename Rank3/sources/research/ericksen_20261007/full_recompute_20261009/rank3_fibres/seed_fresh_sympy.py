"""Export the fresh 17-state system to the original symbolic API names."""
from sage.all import *
from pathlib import Path
import sympy as sy,pickle,json,time
out=Path(__file__).absolute().parent;root=out.parent.parent;st=time.monotonic()
fresh=load(str(out.parent/'foundations/state17_fresh.sobj'))
vs=list(sy.symbols(' '.join(fresh['ring'].variable_names())))
aa,bb=sy.symbols('a b',positive=True);cache={};checks=[]
def cf(c):
 if c not in cache:cache[c]=sy.sympify(str(c).replace('^','**'),locals={'a':aa,'b':bb})
 return cache[c]
def poly(f):return sy.Add(*(cf(c)*sy.Mul(*(v**int(e) for v,e in zip(vs,m) if e)) for m,c in f.dict().items()))
F=[[poly(f) for f in row] for row in fresh['F']]
for i,row in enumerate(F):
 for j,f in enumerate(row):
  reconstructed={tuple(mon):fresh['ring'].base_ring()(str(sy.cancel(co))) for mon,co in sy.Poly(f,*vs).terms() if co != 0}
  assert reconstructed=={tuple(mon):co for mon,co in fresh['F'][i][j].dict().items()},(i,j)
  checks.append('exact coefficient transport '+str((i,j)))
C3=[poly(f) for f in fresh['C3']]
for i,f in enumerate(C3):
 reconstructed={tuple(mon):fresh['ring'].base_ring()(str(sy.cancel(co))) for mon,co in sy.Poly(f,*vs).terms() if co != 0}
 assert reconstructed=={tuple(mon):co for mon,co in fresh['C3'][i].dict().items()},i
 checks.append('exact cubic coefficient transport '+str(i))
api=dict(variables=vs,parameters=(aa,bb),vectorfields=F,constraints_degree3=C3)
target=out/'rebuild/generic_frame_complete_api.pkl';target.parent.mkdir(parents=True,exist_ok=True)
target.write_bytes(pickle.dumps(api))
(out/'fresh_sympy_report.json').write_text(json.dumps(dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st),indent=2))
print('PASS fresh symbolic API',len(checks),time.monotonic()-st,flush=True)
