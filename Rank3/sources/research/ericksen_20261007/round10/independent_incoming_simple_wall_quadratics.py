"""True second coefficients on the two simple normal-kernel walls.

Zero-first-normal source rows require no range correction. Nonzero first
rows are retained as raw input only; no exclusion follows from them alone.
"""
from sage.all import *
from pathlib import Path
import time,resource,os
p=Path(__file__).resolve().parent;st=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(int(Path('/proc/self/statm').read_text().split()[0])*os.sysconf('SC_PAGE_SIZE')+1200*1024**2,)*2)
d=load(str(p/'independent_ratio_incoming12_generic_normal.sobj'));s=load(str(p/'independent_ratio_incoming12_generic_normal_two.sobj'))
U=d['ring'];K=U.base_ring();l,u,v=U.gens();M=d['matrix'];B=s['coefficient_matrix'];exps=s['normal_exponents'];records=[]
print('BASE PAYLOAD',B.nrows(),B.ncols(),'EXPONENTS',len(exps),flush=True)
for name,vec,subs in [('eta',[0,0,0,0,1,0],{u:0}),('k',[-1,0,0,1,0,0],{v:-2})]:
 vals=[]
 weights=[(k,prod(vec[i]**int(e[i]) for i in range(6))) for k,e in enumerate(exps) if sum(e)==2]
 weights=[(k,w) for k,w in weights if w]
 for j in range(M.nrows()):
  coeff=list(B.row(j)) if j<B.nrows() else load(str(p/('independent_ratio_incoming12_generic_normal_two.row%d.sobj'%j)))['coefficients']
  first=U(sum(M[j,k]*vec[k] for k in range(6)).subs(subs));assert not first
  second=U(sum((coeff[k]*w for k,w in weights),U.zero()).subs(subs))
  zero_row=all(not f.subs(subs) for f in M.row(j))
  vals.append({'row':j,'second':second,'zero_first_normal':zero_row})
  if second:
   fac=second.factor();print(name,'row',j,'zero_first',zero_row,'unit',K(fac.unit()).factor(),flush=True)
   print('state',[(str(g),e) for g,e in fac],flush=True)
 records.append({'name':name,'kernel':vec,'values':vals})
save({'ring':U,'items':records,'scope':__doc__},str(p/'independent_incoming_simple_wall_quadratics.sobj'),compress=False)
print('DONE',time.monotonic()-st,flush=True)
