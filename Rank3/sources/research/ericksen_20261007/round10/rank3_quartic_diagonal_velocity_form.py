"""All-parameter six-Q normal form using genuine diagonal ratio derivatives.

Only old h/spectrum/mu guards occur in this exact replacement. No offdiagonal
coordinate is inverted and no state-rank or dimension conclusion is made.
"""
from pathlib import Path
import json,time
from sage.all import load,save
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
d=load(str(out/'joint_fields_direct.sobj'));qs=load(str(base/'round8/rank3_joint_fiber_leading.sobj'))
R=d['ring'];K=R.base_ring();CR=K.ring();a,b,h,z=d['parameters'];x=R.gens();V=d['fields']
A=a*(b-1);C=b*(a-1);D=a-b
k=[-2*(A*z-C*h)/(D*h*z),2*(C-D*z)/(D*z),(A+D*h)/(C*h)]
guard=CR(a*b*(a-1)*(b-1)*D*h*z*(A*z-C*h)*(C-D*z)*(A+D*h))
checks={'old_guards_only':all(not guard.quo_rem(f)[1] for v in k for p in (v.numerator(),v.denominator()) for f,e in p.factor())}
assert checks['old_guards_only']
spec=[(15,0,0,1),(18,0,1,1),(13,1,2,1),(19,1,3,-1),(10,2,4,1),(11,2,5,-1)]
records=[];remainders=[]
for qi,i,xi,sgn in spec:
    Q=R(qs['quartic_source_rows'][qi]);W=V[i,6+i]
    rem=Q-sgn*x[xi]*W/k[i]
    checks['degree_fall_'+str(qi)]=not rem or rem.total_degree()<=3
    assert checks['degree_fall_'+str(qi)]
    remainders.append(rem)
    histogram={}
    for mon in rem.dict():histogram[sum(mon)]=histogram.get(sum(mon),0)+1
    rec={'Q_index':qi,'direction':i,'offdiag':str(x[xi]),'sign':sgn,
         'remainder_degree':int(rem.total_degree()),'remainder_terms':len(rem.dict()),'histogram':histogram}
    records.append(rec);print(json.dumps(rec),flush=True)
save({'ring':R,'parameters':d['parameters'],'diagonal_velocities':[V[i,6+i] for i in range(3)],
      'protected_coefficients':k,'row_specs':spec,'remainders':remainders,'checks':checks,'records':records},str(out/'rank3_quartic_diagonal_velocity_form.sobj'))
(out/'rank3_quartic_diagonal_velocity_form.json').write_text(json.dumps({'status':'PASS_EXACT_ALL_PARAMETER_NORMAL_FORM','coefficient_factors':[str(c.factor()) for c in k],
    'checks':checks,'records':records,'elapsed_seconds':time.monotonic()-start,'scope':__doc__},indent=2)+'\n')
