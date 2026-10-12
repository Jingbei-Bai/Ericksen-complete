"""Exact all-row replay of both endpoint cubic Schur certificates."""
from sage.all import *
from pathlib import Path
import json
OUT=Path(__file__).absolute().parent;checks={};details={}
for tag,pzero in [('q_zero',False),('p_zero',True)]:
 name='homothety_middle_null_endpoint_cubic'+('_pzero' if pzero else '')
 d=load(str(OUT/(name+'.sobj')));S=d['ring'];K=d['parameter_field'];a,r,A2,B2,kfree=K.gens()
 M=d['matrix'];rhs=d['rhs'];C=d['forcing_free_cubic'];rows=d['summary']['rows'];B=M.matrix_from_rows(rows)
 checks[tag+'_full_prior_source_checks']=all(d['source_checks'].values())
 checks[tag+'_matrix_is_whole_76_rows']=M.nrows()==76 and M.ncols()==4 and len(C)==76
 checks[tag+'_all_rows_homogeneous']=all(z==0 for z in rhs)
 expected=-96/r if pzero else -96*r**4
 checks[tag+'_protected_minor']=B.det()==expected
 T=M*B.inverse()
 for i in range(76):
  checks[tag+'_source_row_'+str(i)]=C[i]==sum(T[i,j]*C[rows[j]] for j in range(4))
 checks[tag+'_parameters_unrestricted']=len(K.gens())==5
 assert all(checks.values()),[k for k,v in checks.items() if not v]
 details[tag]={'rows':rows,'determinant':str(expected),'all_four_face_variables':'zero',
  'source_replay_count':sum(d['source_checks'].values())}
out={'status':'PASS_BOTH_ENDPOINT_NEXT_FACES_ONLY','count':len(checks),'checks':checks,'endpoints':details,
 'scope':'Exact cubic source matrices. Does not exclude epsilon=o(delta^2), ramified intermediate orders, or prove no-lift.'}
(OUT/'homothety_middle_null_endpoint_reader.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='checks'}))
