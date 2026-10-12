"""Bridge any early certificate replay to subsequently regenerated sources."""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;root=out.parent.parent;st=time.monotonic();checks=[]
spec={
 'general_curve_extra_quartic_cubics':['rows','leading_rows','diagonal_velocities','diagonal_lower_remainders'],
 'general_curve_remaining_quartics':['rows','leading'],
 'general_curve_homogeneous_normalized':['cubics','quartics','odd_pairs'],
 'general_curve_two_pair_offdiagonal':['cubics','quartics','images'],
 'general_curve_full_amplitude_forms':['A','Q','source_cubics','source_quartics'],
 'general_curve_unit_gradient_directions':['mixed_rows','mixed_differences','curvature_value','old_cubic_rows','old_cubic_scales','amplitude_matrix'],
 'general_curve_centered_six_sources':['sources','linear_matrix'],
 'independent_scalar_lambda_lift':['blocks','source_monodromy','lift_matrix'],
 'independent_scalar_lambda_reduction':['affine_sources','quadratic_source']}
for stem,keys in spec.items():
 fresh=load(str(out/'rebuild/round10'/(stem+'.sobj')));old=load(str(root/'round10'/(stem+'.sobj')))
 for k in keys:
  assert k in fresh,(stem,k,list(fresh))
  assert fresh[k]==old[k],(stem,k)
  checks.append(stem+':'+k)
 print('source equality PASS',stem,flush=True)
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st)
(out/'derived_source_equality.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
