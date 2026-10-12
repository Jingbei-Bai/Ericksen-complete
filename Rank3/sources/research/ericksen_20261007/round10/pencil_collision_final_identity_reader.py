"""Independent exact replay of the explicit three-source protected identity."""
from pathlib import Path
import time,json
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_collision_macaulay_identity.sobj'))
source=load(str(out/'pencil_collision_coefficient_ideal.sobj'))
replay=load(str(out/'pencil_collision_full_source_reader.sobj'))
S=d['guard'].parent();q,T=S.gens();g=q*T*(T-1)*(T-2)*(q-1)*((q-1)*T-4)*(T+4)
checks={}
def check(key,val):checks[key]=bool(val);assert checks[key],key
check('same_original_three_source_factors',list(d['source_polynomials'])==list(source['source_polynomials']))
check('guard_exact',d['guard']==g)
check('power_four',d['power']==4)
check('target_is_guard_power',d['target']==g**4)
check('independent_source_reader_completed',all(replay['checks'].values()))
check('explicit_polynomial_identity',sum((A*P for A,P in zip(d['multipliers'],d['source_polynomials'])),S.zero())==g**4)
check('three_multipliers',len(d['multipliers'])==3)
report={'status':'PASS_INDEPENDENT_EXPLICIT_ORIGINAL_SOURCE_IDENTITY','checks':checks,'check_count':len(checks),
 'multiplier_degrees':[int(v.total_degree()) for v in d['multipliers']],
 'multiplier_terms':[v.number_of_terms() for v in d['multipliers']],
 'elapsed_seconds':time.monotonic()-start,
 'scope':'The source path is now complete: original cubics -> independently replayed t4/t5/t6 -> protected factors -> explicit g^4=sum AiPi. No Groebner-membership oracle is needed to verify this final identity.'}
save({'checks':checks,'identity_difference':sum((A*P for A,P in zip(d['multipliers'],d['source_polynomials'])),S.zero())-g**4},str(out/'pencil_collision_final_identity_reader.sobj'))
(out/'pencil_collision_final_identity_reader.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
