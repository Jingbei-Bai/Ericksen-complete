"""Independent archive replay of the all-parameter physical infinity bridge.

This is NOT a no-common-root certificate. It validates the three necessary and
sufficient compatibility conics after the protected third-jet elimination.
"""
from sage.all import *
from pathlib import Path
import time,signal,json,resource
O=Path(__file__).resolve().parent;st=time.time();signal.alarm(150)
resource.setrlimit(resource.RLIMIT_CORE,(0,0))
vm=int(Path('/proc/self/statm').read_text().split()[0])*4096
resource.setrlimit(resource.RLIMIT_AS,(vm+1300*1024**2,vm+1300*1024**2))
checks={}
sym=load(str(O/'rr_cartesian_twelve_symbol_certificate.sobj'))
for tag in ['plus','minus']:
 d=load(str(O/f'rr_null_physical_protected_graph_{tag}.sobj'))
 old=load(str(O/f'rr_null_physical_hessian_source_{tag}.sobj'))
 src=load(str(O/f'rr_null_physical_late_full_source_{tag}.sobj'))
 raw=load(str(O/f'rr_null_physical_before_graph_quadratics_{tag}.sobj'))
 out=load(str(O/f'rr_null_physical_quadratics_{tag}.sobj'))
 K=d['field'];P=K.ring();T=d['ring'];aa,bb,chi=T.gens();s,t,p,q,c=K.gens();J=s*t*c
 h=old['ring'].hom([aa,bb]+d['graph'],T)
 for n,f in enumerate(old['first_equations']):checks[f'{tag}_first_source_{n}']=h(f)==0
 checks[f'{tag}_source_axes_null']=d['axes']*d['axes']==0
 checks[f'{tag}_target_axes_null']=d['axet']*d['axet']==0
 baseM=sym['matrix'].change_ring(P)
 scaled=diagonal_matrix(P,[J if i in [2,5,8] else 1 for i in range(12)])*baseM
 checks[f'{tag}_protected_twelve_symbol']=src['matrix'][:12,:]==scaled
 ker=raw['kernel']
 checks[f'{tag}_all_fifteen_symbol']=ker*src['matrix']==0
 checks[f'{tag}_last_three_identity']=ker[:,12:]==identity_matrix(K,3)
 checks[f'{tag}_raw_all_source_contraction']=list(ker*src['offset'])==raw['quadratics_raw_H']
 ph=raw['ring'].hom([aa,bb]+d['graph'],T)
 for n,f in enumerate(raw['quadratics_raw_H']):
  checks[f'{tag}_conic_{n}_source_replay']=ph(f)==out['raw_quadratics'][n]
 checks[f'{tag}_no_removed_content']=out['raw_quadratics']==out['quadratics'] and out['multipliers']==[K(1)]*3
 checks[f'{tag}_all_amplitudes_retained']=all(f.total_degree()==2 and f.is_homogeneous() for f in out['quadratics'])
 print(tag,sum(checks.values()),len(checks),time.time()-st,flush=True)
assert all(checks.values()),checks
report={'status':'PASS_COMPLETE_SOURCE_BRIDGE_ONLY_NO_PROJECTIVE_EXCLUSION',
 'checks':checks,'count':len(checks),'seconds':time.time()-st,
 'domain':'s,t,c>0; p*q*(s^2-t^2)!=0; alpha,beta,chi arbitrary complex, including zero',
 'note':'The original source pressures and angular norms are fixed finite data; their highest homogeneous parts are zero. The conics still require a whole-parameter no-common-root proof.'}
(O/'rr_null_physical_source_bridge_reader.json').write_text(json.dumps(report,indent=2),encoding='utf8')
print(json.dumps(report,indent=2),flush=True)
