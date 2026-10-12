from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;old=out.parent.parent/'round10';st=time.monotonic();checks=[]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
entries=[
('independent_delta_critical_canonical_amplitude','independent_delta_critical_canonical_amplitude',['radical_squares','extension_matrix','row_phases','column_phases','material_matrix','old4_determinant']),
('independent_supercritical_exceptional_even','independent_supercritical_exceptional_even',['even','root_e_phases']),
('independent_supercritical_exceptional_even','independent_supercritical_line_even',['normal_matrix','determinant','e_graph','even_pullbacks']),
('independent_delta_subcritical_H_projected_normals','independent_delta_subcritical_H_projected_normals',['parameters','p21_value','e_graph','C','C_determinant','C_adjugate','R_p23','R_kappa','R_delta','Q0_e_row','Q0_p23','Q0_kappa','Q0_delta','Q0_line_schur','Q0_kappa_schur','Q0_delta_schur','items']),
('independent_supercritical_normal_direction','independent_supercritical_normal_direction',['sources','items','groebner_basis','target_lifts']),
]
for creator,name,keys in entries:
 path=out/'source_replays'/creator/'generated/round10'/(name+'.sobj')
 a=load(str(path));b=load(str(old/(name+'.sobj')))
 for key in keys:ck(name+' full '+key,a[key]==b[key])
 if name=='independent_supercritical_normal_direction':
  for target,coeffs in a['target_lifts']:ck('fresh literal coefficient-ideal witness',sum(c*f for c,f in zip(coeffs,a['sources']))==target)
 print('BRIDGED',name,flush=True)
(out/'fresh_canonical_projection_links.json').write_text(json.dumps(dict(status='PASS',count=len(checks),checks=checks,seconds=time.monotonic()-st,source='Current isolated producers reconstruct every displayed source array from already fresh-linked critical22x5/all3even and actual two-parameter jets; all downstream consumed arrays compared exactly.',producers=[dict(script=a,object=b,checked_keys=c) for a,b,c in entries]),indent=2))
print('PASS',len(checks),time.monotonic()-st,flush=True)
