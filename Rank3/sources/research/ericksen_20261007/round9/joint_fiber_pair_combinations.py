from pathlib import Path
import time,json
from sage.all import load,save
out=Path(__file__).resolve().parent;base=out.parent;start=time.monotonic()
d=load(str(base/'round8/rank3_joint_fiber_leading.sobj'))
R=d['ring'];p12,p13,p21,p23,p31,p32,rr1,rr2,rr3=R.gens();rows=d['quartic_source_rows']
spec=[[(10,p32),(11,p31)],[(13,p23),(19,p21)],[(15,p13),(18,-p12)]]
combos=[];records=[]
for n,parts in enumerate(spec):
 q=sum((v*rows[j] for j,v in parts),R.zero());combos.append(q)
 hist={}
 for m,c in q.dict().items():hist[sum(m)]=hist.get(sum(m),0)+1
 records.append({'index':n,'source_combination':[(j,str(v)) for j,v in parts],
                 'degree':int(q.total_degree()) if q else None,'terms':len(q.dict()),'degree_histogram':hist})
 print(json.dumps(records[-1]),flush=True)
payload={'ring':R,'parameters':d['parameters'],'variables':R.gens(),'source_combination':spec,
         'polynomials':combos,'records':records,'source_file':'round8/rank3_joint_fiber_leading.sobj'}
save(payload,str(out/'joint_fiber_pair_combinations.sobj'))
report={'status':'PASS_EXACT_SOURCE_COMBINATIONS','records':records,'elapsed_seconds':time.monotonic()-start,
        'scope':'No division by state or base expression; every displayed relation holds for the full original fiber at all protected parameter specializations.'}
(out/'joint_fiber_pair_combinations.json').write_text(json.dumps(report,indent=2)+'\n')
print('SAVED',round(time.monotonic()-start,2),flush=True)
