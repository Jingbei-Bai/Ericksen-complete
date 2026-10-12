"""Eliminate N from the small lambda lift using an inherited guard only."""
from pathlib import Path
import time,json
from sage.all import load,save,QQ,PolynomialRing,matrix,vector

out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'independent_scalar_lambda_lift.sobj'))
O=d['ring'];oldH,oldZ,oldS,oldN=O.gens()
R=PolynomialRing(QQ,names=('H','Z','S','lambda','tau'));H,Z,S,lam,tau=R.gens()
def embed(c):return R(c(H,Z,S,0))
L=matrix(R,4,4,[embed(c) for c in d['lift_matrix'].list()])
e=L[3,2];c0=L[3,0];d0=-L[3,1];d1=-L[3,3]
checks={'protected_N_pivot':e==-2*H*Z*(H-Z)}
Nnum=d0*lam+d1*tau-c0;Nvalue=R.fraction_field()(Nnum/e)
affine=[L[i,0]+L[i,1]*lam+L[i,3]*tau for i in range(3)]
quad=e*tau-d0*lam**2-d1*lam*tau+c0*lam
checks['product_constraint']=R.fraction_field()(quad)==e*(tau-lam*Nvalue)
checks['fourth_row_zero']=c0-d0*lam+e*Nvalue-d1*tau==0
checks['first_rows_independent_N']=all(L[i,2]==0 for i in range(3))
K=R.fraction_field()
actual_eq=[K(d['equations'][i].list()[0](H,Z,S,Nvalue))+
           lam*K(d['equations'][i].list()[1](H,Z,S,Nvalue)) for i in range(4)]
# Actual equations differ from lifted equations only by the controlled
# product defect. This proves both directions in the localization at e.
for i in range(3):
    checks['original_equation_mod_product_'+str(i)]=actual_eq[i]-affine[i]==L[i,3]*(lam*Nvalue-tau)
checks['original_fourth_mod_product']=actual_eq[3]==L[3,3]*(lam*Nvalue-tau)

T=matrix(R,3,3,[[L[i,0],L[i,1],L[i,3]] for i in range(3)])
compat=T.det()
checks['affine_matrix_replay']=T*vector(R,[1,lam,tau])==vector(R,affine)
assert all(checks.values()),[k for k,v in checks.items() if not v]
payload={'ring':R,'variables':R.gens(),'affine_sources':affine,'quadratic_source':quad,
 'N_numerator':Nnum,'N_denominator':e,'N_reconstruction':Nvalue,'affine_matrix':T,
 'affine_compatibility':compat,'affine_compatibility_factorization':compat.factor(),
 'checks':checks,'source_file':'independent_scalar_lambda_lift.sobj',
 'scope':'Exact equivalence to the four lambda-lift equations after adjoining tau=lambda*N, with N eliminated using only the old H*Z*(H-Z) guard. Lambda=0 and D00=0 are fully retained. Scalar-monodromy equivalence still needs det(B01)*det(B02)!=0.'}
save(payload,str(out/'independent_scalar_lambda_reduction.sobj'))
report={'status':'EXACT_PROTECTED_REDUCTION_PASS','checks':checks,
 'affine_terms':[len(p.dict()) for p in affine],'quadratic_terms':len(quad.dict()),
 'compatibility_terms':len(compat.dict()),'compatibility_factorization':str(compat.factor()),
 'elapsed_seconds':time.monotonic()-start,'scope':payload['scope']}
(out/'independent_scalar_lambda_reduction.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
