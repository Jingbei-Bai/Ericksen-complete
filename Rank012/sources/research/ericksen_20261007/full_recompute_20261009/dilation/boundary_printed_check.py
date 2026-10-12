"""Check every printed shared Section 7/8 normal pivot against regenerated sources."""
from sage.all import *
from pathlib import Path
import json,time
HERE=Path(__file__).absolute().parent;D=HERE/'rebuild'/'round12';started=time.time();checks={};records={}
onezero=load(str(D/'homothety_onezero_fiber_leading.sobj'))
Ka=onezero['ring'].base_ring();aa,bb,vv=Ka.gens()
for i,rec0 in enumerate(onezero['records']):
 expected_obs=((aa-bb)**3*vv**4/(4*aa**3*bb**3) if i==0 else (bb-aa)*(bb-1)**2*vv**4/(4*aa*bb**3))
 expected_T=(2*(aa-bb)*vv**4/(aa*(bb-1)) if i==0 else -2*vv**4)
 expected_lead=(-vv*bb*(aa-1)/(aa-bb) if i==0 else vv*aa*bb*(bb-1)*(aa-1)/(aa-bb)**2)
 checks['single_zero_observation_printed_'+str(i)]=rec0['observation_matrix'].matrix_from_columns(list(range(1,9))).det()==expected_obs
 checks['single_zero_T_printed_'+str(i)]=rec0['t_pivot']==expected_T
 checks['single_zero_leader_printed_'+str(i)]=rec0['leaders'][1]==expected_lead
axis=load(str(D/'homothety_axis_jet_bound_reader.sobj'));Ka=axis['ring'].base_ring();aa,bb,ss,tt=Ka.gens()
checks['axis_T_printed_minor']=axis['t_matrix'].matrix_from_rows_and_columns(axis['t_rows'],axis['t_columns']).det().subs({ss:0,tt:0})==(aa-bb)/(aa*(bb-1))
checks['axis_U_identity_block']=axis['u_matrix'].matrix_from_rows(axis['u_rows'])==identity_matrix(Ka,2)
checks['axis_quartic_index_is_Q6']=axis['Q4_index']==6
normal=load(str(D/'homothety_axis_plane_normal_sources.sobj'))
forcing=load(str(D/'homothety_axis_plane_forcing_schur.sobj'))
quad=load(str(D/'homothety_axis_plane_raw_quadratic_schur.sobj'))
K=normal['parameter_field'];a,b,p,q=K.gens()
L=forcing['forcing_matrix'];Q=quad['quadratic_matrix'];rangeM=normal['observation_plus_normal_matrix']
wrong_rows=list(range(6))+[8,10,13,14,15,17,19];cc=[0,1,3,4,6,7,10,11,12,13,14,15,16]
checks['dil:original_six_angular_selection_is_singular']=rangeM.matrix_from_rows_and_columns(wrong_rows,cc).det()==0
rr=[0,1,2,5,6,7,8,10,13,14,15,17,19]
checks['dil:corrected_range_minor']=rangeM.matrix_from_rows_and_columns(rr,cc).det()==(a-1)*(a-b)/(4*a*a*(b-1))
records['range_selection_correction']={'original_observations':[0,1,2,3,4,5],
    'correct_observations':[0,1,2,5,6,7],
    'meaning':'A1,A2,A3,B3,K1,K2 (four angular and two pressure observations)',
    'homothety_rows':[0,2,5,6,7,9,11], 'columns':cc,
    'original_determinant':'0','corrected_determinant':'(a-1)*(a-b)/(4*a^2*(b-1))'}
checks['iso:scale_column_in_pressure_image']=L.column(8)==L.column(6)+L.column(7)
forcing_target=matrix.zero(K,5,9)
for j,k in enumerate([0,1,3,4]):forcing_target[j,k]=1
forcing_target[4]=vector(K,[0,0,-2*a*p*q,0,0,2*a*b*p*q,-a*(b-1)*p*p,-(a-b)*q*q,-a*(b-1)*p*p-(a-b)*q*q])
checks['dil:entire_forcing_kernel']=L.row_space()==forcing_target.row_space()
def projected(L0,Q0,rr,cc):
 C=L0.matrix_from_columns(cc)*L0.matrix_from_rows_and_columns(rr,cc).inverse()
 assert C*L0.matrix_from_rows(rr)==L0
 return Q0-C*Q0.matrix_from_rows(rr)
I=[0,2,5,7,3,8,9]
rowsets=[[12,16,18,20,34,46,73],[11,12,18,21,34,47,75],[12,16,18,21,34,46,73],[12,16,18,22,34,46,73]]
for mode in ['general','qzero','wall','doublewall','pzero']:
 if mode=='pzero':
  sub={p:0};frows=[3,4,11,21,27];fcols=[0,1,3,4,7]
 else:
  sub={p:1};frows=[3,4,9,11,27];fcols=[0,1,2,3,4]
  if mode=='qzero':sub[q]=0;frows=[3,4,9,22,27];fcols=[0,1,3,4,6]
 QQ0=Q.apply_map(lambda x:x.subs(sub));LL=L.apply_map(lambda x:x.subs(sub))
 PP=projected(LL,QQ0,frows,fcols)
 if mode in ['wall','doublewall']:
  asub=b*q*q/(q*q+1-b)
  PP=PP.apply_map(lambda x:x.subs({a:asub}))
  if mode=='doublewall':PP=PP.apply_map(lambda x:x.subs({b:(7*q*q+3)/4}))
 if mode=='general':
  rows=rowsets[0];expected=-49152*(b-a)*(b-1)**2*a*a*q**5*((b-a)*q*q+a*(b-1))/(b**4*(a-1)**4)
 elif mode=='qzero':rows=rowsets[1];expected=41472*(b-a)*(b-1)**4*a**4/(b**5*(a-1)**5)
 elif mode=='wall':rows=rowsets[2];expected=2048*q**10*(7*q*q-4*b+3)/(b*(q*q+1)**3)
 elif mode=='doublewall':rows=rowsets[3];expected=-49152*q**9/((q*q+1)**2*(3*q*q-1))
 else:
  # Select a sparse row basis over the exact parameter field, then verify
  # the endpoint determinant printed without a row selector in the paper.
  rows=[];rank=0
  order=sorted(range(PP.nrows()),key=lambda i:sum(len(x.numerator().dict())+len(x.denominator().dict()) for x in PP.row(i) if x))
  for i in order:
   newrank=PP.matrix_from_rows(rows+[i]).rank()
   if newrank>rank:rows.append(i);rank=newrank
  expected=4608*(a-b)**5*q**10/(b**4*(a-1)**4*(b-1)*a)
 cols=I if mode!='pzero' else list(PP.pivots())
 block=PP.matrix_from_rows_and_columns(rows,cols);det=block.det()
 if expected is not None:checks['dil:quadratic_printed_'+mode]=det==expected
 checks['dil:quadratic_complete_rank7_'+mode]=len(cols)==7 and PP==PP.matrix_from_columns(cols)*block.inverse()*PP.matrix_from_rows(rows)
 records['quadratic_'+mode]={'rows':rows,'columns':cols,'minor':str(det.factor())}
 if mode=='general':
  target=matrix.zero(K,7,10)
  target[0,0]=1;target[0,4]=(b-a)/(a*(b-1))
  target[1,2]=1;target[2,5]=1;target[3,7]=1
  target[4,3]=1;target[4,6]=(a-b)*q/(a*(b-1))
  target[5,8]=1;target[6,9]=1
  checks['dil:printed_seven_generator_lift']=target.matrix_from_columns(I)*block.inverse()*PP.matrix_from_rows(rows)==target

cubic=load(str(D/'homothety_middle_null_cubic_solution.sobj'))
K1=cubic['matrix'].base_ring();aa,pp,qq,r,A2,B2,k1=K1.gens()
rows=[46,20,48,12];M=cubic['matrix'];rhs=cubic['rhs']
checks['dil:middle_four_printed_minor']=M.matrix_from_rows(rows).det()==1024*r*qq**3*pp**3*(pp*r+qq)**3/(1+r*r)
center=vector(K1,[0,0,0,-2*r*pp*pp*qq*(1+r*r)/(pp*r+qq)])
checks['dil:middle_moving_cubic_center']=M*center==rhs
for j in [1,2,3]:checks['dil:middle_observation_free_column_'+str(j)]=all(x.derivative(v)==0 for x in M.column(j) for v in [A2,B2,k1])
quartic=load(str(D/'homothety_middle_null_quartic_solution.sobj'))
M4=quartic['matrix'];rhs4=quartic['rhs']
sol=M4.matrix_from_rows(rows).inverse()*vector(K1,[rhs4[i] for i in rows]);res=M4*sol-rhs4
checks['dil:middle_quartic_first_obstruction']=res[16]==8*pp*qq*(pp*r-qq)*(1+r*r)**2/(pp*r+qq)
checks['dil:middle_quartic_last_obstruction']=res[34].subs({qq:pp*r})==8*r*r*pp**4*(1+r*r)**2
for pzero in [False,True]:
 suffix='_pzero' if pzero else '';data=load(str(D/('homothety_middle_endpoint_uniform_quartic'+suffix+'.sobj')))
 K2=data['field'];aa,r,A2,B2,kfree,v=K2.gens();rr=[9,12,16,46] if pzero else [11,12,16,47]
 FP=data['forcing'].matrix_from_rows_and_columns(data['forcing_rows'],data['forcing_columns'])
 expected=(r*r+aa)**2/(aa*aa*r*r) if pzero else (r*r+aa)**2/r**4
 checks['dil:endpoint_forcing_minor'+suffix]=FP.det().subs({v:0})==expected
 expected=-96/r if pzero else -96*r**4
 checks['dil:endpoint_four_minor'+suffix]=data['matrix'].matrix_from_rows(rr).det().subs({v:0})==expected
 expected=20*(1+r*r) if pzero else 20*r*r*(1+r*r)
 checks['dil:endpoint_quartic_obstruction'+suffix]=data['residuals'][22 if pzero else 21].subs({v:0})==expected

resonance=load(str(D/'homothety_middle_resonant_whole_coefficients.sobj'))
checks['dil:resonance_every_coefficient']=all(resonance['checks'].values())
rec={'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'count':len(checks),'records':records,'seconds':time.time()-started}
(HERE/'boundary_printed_check.json').write_text(json.dumps(rec,indent=2))
print(json.dumps({'status':rec['status'],'count':len(checks),'failed':[n for n,v in checks.items() if not v],'seconds':rec['seconds']}),flush=True)
assert all(checks.values()),[n for n,v in checks.items() if not v]
