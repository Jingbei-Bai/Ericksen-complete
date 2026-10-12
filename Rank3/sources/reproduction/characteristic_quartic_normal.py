"""Invert the diagonal u coefficients of six complete quartic sources."""
from pathlib import Path
import argparse,json,pickle,time
from sage.all import load,matrix,identity_matrix,lcm

parser=argparse.ArgumentParser()
parser.add_argument('--work-root',type=Path,required=True)
args=parser.parse_args();root=args.work_root.resolve();started=time.monotonic()
fresh=load(str(root/'research/ericksen_20261007/full_recompute_20261009/foundations/state17_fresh.sobj'))
P=fresh['ring'];K=P.base_ring();a,b=K.gens();x=P.gens();Q=fresh['Q4'];checks=[]
def ck(name,value):
    assert value,name
    checks.append(name)
indices=[0,1,3,5,6,12]
targets=[x[j]*x[15+i] for j in (0,4,8) for i in (0,1)]
ck('complete quartic source and state dimensions',P.ngens()==17 and len(Q)==21)
selected=[Q[j] for j in indices]
B=matrix(K,6,6,lambda i,j:selected[i].monomial_coefficient(targets[j]))
det=2464*(a-b)*b**2*(a-1)**2/(a**3*(b-1)**3)
ck('complete determinant has only admissible spectral factors',B.det()==det)
L=B.inverse()
ck('exact left and right inverse',L*B==B*L==identity_matrix(K,6))
common_denominator=lcm([c.denominator() for c in L.list()])
for factor,exponent in common_denominator.factor():
    if factor.is_constant():continue
    ratios=[K(factor)/guard for guard in (a,b,a-1,b-1,a-b)]
    ck('inverse denominator factor '+str(factor),any(c.numerator().is_constant() and c.denominator().is_constant() for c in ratios))
for i in range(6):
    combined=sum((L[i,j]*selected[j] for j in range(6)),P.zero())
    u_part=P({tuple(e):c for e,c in combined.dict().items() if any(e[15:])})
    ck('complete quartic u part '+str(i),u_part==targets[i])
payload={'coefficient_transform':L,'source_quartic_indices':indices}
target=root/'research/ericksen_20261007/round8/rank3_k_quartic_diagonal_normal.pkl'
target.parent.mkdir(parents=True,exist_ok=True)
target.write_bytes(pickle.dumps(payload))
ck('exact interface readback',pickle.loads(target.read_bytes())==payload)
report={'status':'PASS','count':len(checks),'checks':checks,'seconds':time.monotonic()-started}
(target.parent/'fresh_characteristic_quartic_normal.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS full diagonal quartic normal form',len(checks),round(report['seconds'],2),flush=True)
