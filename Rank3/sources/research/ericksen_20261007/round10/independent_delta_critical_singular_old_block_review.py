"""Audit the singular old 4-block on all six nontrivial direction walls."""
from sage.all import *
from pathlib import Path
import time,json
p=Path(__file__).resolve().parent;st=time.monotonic();d=load(str(p/'independent_critical_exceptional_old_block.sobj'));K=d['field'];a,r,H,Z,B,C=K.gens();q=1+r;out=[]
for it in d['records'][1:]:
 s=K(it['p21']);target=(B*B+q*s*s) if it['direction_type']=='omega_omega' else (r*C*C+q)
 quotient=it['restricted_block']/target;assert quotient*target==it['restricted_block'];factors=quotient.factor()
 print('CASE',it['direction_type'],'FACTORS',[(str(f),int(e)) for f,e in factors],flush=True)
 out.append({'direction_type':it['direction_type'],'canonical_target':target,'quotient':quotient,'factorization':factors})
save({'source_file':'independent_critical_exceptional_old_block.sobj','field':K,'number_field':d['number_field'],'records':out,'scope':__doc__},str(p/'independent_delta_critical_singular_old_block_review.sobj'),compress=False)
(p/'independent_delta_critical_singular_old_block_review.json').write_text(json.dumps({'cases':len(out),'seconds':time.monotonic()-st,'scope':__doc__},indent=2)+'\n')
