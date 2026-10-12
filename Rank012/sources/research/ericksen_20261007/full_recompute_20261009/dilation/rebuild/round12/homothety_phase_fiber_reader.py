from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).absolute().parent;d=load(str(OUT/'homothety_invariant_linear_chart.sobj'));o=d['source'];R=o['ring'];K=R.base_ring();a,b,s,t=K.gens();M=d['matrix'];v=d['fiber_kernel'];ld=d['quadratic_fiber_leaders'];cc=load(str(OUT/'homothety_fiber_leading_compression.sobj'));ABC=cc['F'].parent();X,Y,Z=ABC.gens();F=cc['F'];G=cc['G'];checks={}
checks['kernel_source']=M*v==vector(K,[0]*8)
checks['kernel_q_normalized']=v[2]==1
checks['first_t_protected_pivot']=o['determinant']==2*(a-b)*s*s*t/(a*(b-1))
d2=load(str(OUT/'homothety_second_jet_reduction.sobj'))
checks['highest_u_protected_pivot']=d2['u_pivot_det']==-b*(a-1)/(2*a*(b-1))
Qm=b*(a-1)*s*s*t*t+a*(b-1)*s*s-(a-b)*t*t
idx=[i for i in range(9) if i!=2];BR=K.ring();N=matrix(BR,8,8,[BR(2*a*b*z) for z in M.matrix_from_columns(idx).list()]);det=K(N.det())/(2*a*b)**8
checks['observation_pivot']=det==(b-a)*s*(b-1)*Qm/(4*a*a*b**3)
checks['quadratic_source_leaders_proportional']=ld[2]==b*(a-1)*ld[0]/(a*(b-1))
xx=b*(a-1)*s*s*t*t;yy=a*(b-1)*s*s;zz=(a-b)*t*t
ff=K(F(xx,yy,zz));gg=K(G(xx,yy,zz))
fac0=-(b*(a-1)*s*s+a-b)/((b-a)**3*s**3*t*(b-1)*a*Qm)
fac1=b*(a-1)*(b*(a-1)*t*t+a*(b-1))/((b-a)**2*(b-1)**2*a*a*s*Qm**2)
checks['first_full_leader_compact']=ld[0]==fac0*gg
checks['second_full_leader_compact']=ld[1]==fac1*ff*ff
checks['compact_difference']=G-F*F==-8*(-X+Y-Z)*Z*(Y-Z)*(X+Y+Z)
checks['first_prefactor_wall_excluded']=F.subs({X:-Z})==Y*Y
checks['G_factor_wall_one']=F.subs({Y:X+Z})==4*X*(X+Z)
checks['G_factor_wall_two']=F.subs({Z:-X-Y})==-4*Y*(X+Y)
assert all(checks.values()),checks
out={'status':'PASS','checks':checks,'check_count':len(checks),'scope':'Entire real nonzero-component homothety chart, all positive distinct material spectra and all real source/target twists and longitudinal pressures; no assertion of full pairing exclusion or global uniform bound.'}
(OUT/'homothety_phase_fiber_reader.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
