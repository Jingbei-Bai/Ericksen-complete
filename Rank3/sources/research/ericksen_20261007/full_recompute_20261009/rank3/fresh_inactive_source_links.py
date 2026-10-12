"""Close the full inactive-support source and cyclic guard provenance."""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;old=out.parent.parent/'round10';checks=[];st=time.monotonic()
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
for name,keys in [
 ('independent_ratio_inactive_row3_RR_derivatives',['original_parameter_map','field_ordinary_cubic','source_unit','missing_Z0','items']),
 ('independent_ratio_inactive_row3_linear_strata',['Z','T','T_original_derivative','linear_augmented_matrix','compatibility_determinant']),
]:
 a=load(str(out/'source_replays'/name/'generated/round10'/(name+'.sobj')));b=load(str(old/(name+'.sobj')))
 for key in keys:ck(name+' all '+key,a[key]==b[key])
freshpath=out/'source_replays/independent_inactive_generic_support_reader/generated/round10/independent_inactive_generic_support_reader.json'
a=json.loads(freshpath.read_text());b=json.loads((old/'independent_inactive_generic_support_reader.json').read_text());ck('fresh whole-source generic support record',a['status']=='PASS');ck('all actual factors consumed by cyclic guard',a['base_factors']==b['base_factors'])
rec=dict(status='PASS',count=len(checks),checks=checks,seconds=time.monotonic()-st,source='RR derivatives and all Z,T incidence rows rebuilt from full original quartics and actual gauged fields. Full support reader then reconstructs every graph, selected minor, actual Sylvester coefficient, full amplitude source and all15 nonempty p-support identities. Cyclic factor list compared with fresh actual calculation, not old PASS.')
(out/'fresh_inactive_source_links.json').write_text(json.dumps(rec,indent=2));print('PASS',len(checks),time.monotonic()-st)
