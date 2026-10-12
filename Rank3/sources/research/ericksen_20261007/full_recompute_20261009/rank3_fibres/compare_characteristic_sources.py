"""Complete equality bridge for every early characteristic witness input."""
from sage.all import *
from pathlib import Path
import json,time
out=Path(__file__).absolute().parent;root=out.parent.parent;st=time.monotonic();checks=[]
spec={
 'round9/rank1_ratio_characteristic_quartics_partial':['t_numerators_over_p23','Q_scaled_lower','Q_u_coefficients'],
 'round9/rank1_ratio_characteristic_second_tangency':['scaled_lower_rows','u_coefficient_matrix'],
 'round9/rank1_ratio_characteristic_second_graph':['u_numerators_over_p23','remaining_second_tangencies','selected_quartic_rows'],
 'round9/rank1_ratio_characteristic_linear_pair':['linear_source_rows','remaining_mixed_source','quartic_source_Q0'],
 'round10/pencil_fast_numerators':['rows','determinant','cramer_numerators'],
 'round10/pencil_vzero_sone_cubics':['cubics','c6','d62','d60'],
 'round10/pencil_aqminusone_newton':['cubics'],
 'round10/pencil_aqminusone_upper_walls':['coefficients','remainders'],
 'round10/pencil_aqminusone_wall_factors':['N','J'],
 'round10/pencil_aqminusone_Azero':['coefficient_Hb109','remainder'],
 'round10/pencil_aqminusone_Ttwo_sylvester':['coefficients']}
for stem,keys in spec.items():
 fresh=load(str(out/'rebuild'/(stem+'.sobj')));old=load(str(root/(stem+'.sobj')))
 for k in keys:
  assert k in fresh,(stem,k,list(fresh))
  assert fresh[k]==old[k],(stem,k)
  checks.append(stem+':'+k)
 print('source equality PASS',stem,flush=True)
report=dict(status='PASS',checks=checks,count=len(checks),seconds=time.monotonic()-st)
(out/'characteristic_source_equality.json').write_text(json.dumps(report,indent=2));print(report,flush=True)
