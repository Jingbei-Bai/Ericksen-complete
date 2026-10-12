"""Verify the exact kernel and minors PRINTED in Section 8 against fresh jets."""
from sage.all import *
from pathlib import Path
import json,time
HERE=Path(__file__).absolute().parent; OUT=HERE/'rebuild'/'round23';started=time.time()
checks={};records={}
for mode in ['general','pzero','qzero','resonance']:
 data=load(str(OUT/('dilation_state_column_linear_module_'+mode+'.sobj')))
 K=data['field'];gg=K.gens();a=gg[0]
 if mode in ['general','resonance']:
  zz=gg[-1];p=2*zz/(1+zz*zz);q=(1-zz*zz)/(1+zz*zz)
  b=gg[1] if mode=='general' else a*p*p+q*q
 else:b=gg[1];p=K(mode=='qzero');q=K(mode=='pzero')
 L=data['ambient_jacobians'];v=vector(K,[q,-p,0]);ax=vector(K,[0,0,1]);w=ax.cross_product(v)
 DB=sum((ax[i]*L[i] for i in range(3)),matrix.zero(K,23))
 raw=DB.stack(data['raw_source_jacobian'])
 if mode=='pzero':
  rr=[0,2,3,5,6,7,9,11,13,18,25,74,75,76]
  cc=[0,3,6,7,10,11,12,13,14,15,16,17,19,22]
  free=[1,2,4,5,8,9,18,20,21]
  expected=2*(a-1)**2/(a*a*(b-1)**2)
 elif mode=='qzero':
  rr=[0,1,2,5,6,7,9,11,14,17,25,74,75,76]
  cc=[0,1,6,7,10,11,12,13,14,15,16,18,19,22]
  free=[2,3,4,5,8,9,17,20,21]
  expected=2*(b-a)*(a-1)**2/(a*(b-1)**3)
 else:
  rr=[0,1,2,3,5,6,7,9,11,17,25,74,75,76]
  cc=[0,1,6,7,10,11,12,13,14,15,16,17,19,22]
  free=[2,3,4,5,8,9,18,20,21]
  expected=2*p**4*q*q*(a-1)**2/(a*(b-1)**2)
 checks[mode+'_printed_minor']=raw.matrix_from_rows_and_columns(rr,cc).det()==expected
 N=matrix.zero(K,23,9)
 for k,j in enumerate(free):N[j,k]=1
 N[0,2]=-1;N[0,4]=1;N[0,7]=2*p;N[0,8]=2*q
 if mode=='pzero':N[3,7]=1;N[12,7]=1;N[19,7]=-1
 elif mode=='qzero':N[1,8]=1;N[10,8]=-1;N[19,8]=1
 else:
  g=b*(a-1)/(a*(b-1));d=(a-b)/(a*(b-1))
  N[1,1]=d*q*q/(p*p);N[1,2]=-g*q/p;N[1,4]=q/p
  N[1,7]=q*(1-g*q*q)/(p*p);N[1,8]=(1+g*q*q)/p
  N[10,7]=p*q;N[10,8]=-p*p
  N[12,7]=q*q;N[12,8]=-p*q
  N[15,7]=2*p*q*q;N[15,8]=-2*p*p*q
  N[17,6]=p/q;N[19,7]=-q;N[19,8]=p
 checks[mode+'_printed_kernel']=raw*N==0
 for i in range(3):checks[mode+'_complete_linear_closure_'+str(i)]=raw*L[i]*N==0
 R=(sum((v[i]*L[i] for i in range(3)),matrix.zero(K,23))*N).matrix_from_rows(free)
 A=(sum((w[i]*L[i] for i in range(3)),matrix.zero(K,23))*N).matrix_from_rows(free)
 Id=identity_matrix(K,9)
 P=Id+5*R/3-R**2/12-7*R**3/6-5*R**4/12
 checks[mode+'_printed_characteristic']=R.charpoly()==R.charpoly().parent()('x*(x+2)*(x-1)^2*(x+1)^5')
 checks[mode+'_printed_annihilator']=R*(R-Id)*(R+2*Id)*(R+Id)**2==0
 checks[mode+'_commute']=R*A==A*R
 checks[mode+'_projector']=P*P==P
 checks[mode+'_projector_trace']=P.trace()==3
 checks[mode+'_angular_vanish']=A*P==0
 checks[mode+'_semisimple']=R*(R-Id)*P==0
 records[mode]={'minor':str(expected.factor()),'R':[[str(x) for x in row] for row in R.rows()],
                'A':[[str(x) for x in row] for row in A.rows()], 'free_rows':free}
assert all(checks.values()),[n for n,v in checks.items() if not v]
rec={'status':'PASS_PRINTED_KERNEL_MINORS_AND_SPECTRA','checks':checks,'count':len(checks),
     'records':records,'seconds':time.time()-started,
     'source':'Freshly differentiated original fields; no saved carrier matrices loaded.'}
(HERE/'carrier_printed_check.json').write_text(json.dumps(rec,indent=2))
print(json.dumps({'status':rec['status'],'count':len(checks),'seconds':rec['seconds']}),flush=True)
