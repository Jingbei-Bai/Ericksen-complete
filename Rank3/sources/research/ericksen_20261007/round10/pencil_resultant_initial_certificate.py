"""First exact resultant coefficients, by truncated Sylvester determinant."""
from pathlib import Path
import json,time
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_newton_common_root.sobj'));F=d['series'][0];G=d['series'][1]
K=F[0].base_ring();C=K.ring();a,q,T=C.gens();X=F[0].parent().gen();r=q-1
vecF=[[C(F[k][j]) for k in range(3)] for j in reversed(range(4))]
vecG=[[C(G[k][j]) for k in range(3)] for j in reversed(range(4))]
zero=[C.zero()]*3
M=[]
for vv in (vecF,vecG):
    for shift in range(3):M.append([zero]*shift+vv+[zero]*(2-shift))
def mul(x,y):return [sum((x[j]*y[k-j] for j in range(k+1)),C.zero()) for k in range(3)]
dp={0:[C.one(),C.zero(),C.zero()]}
for row in range(6):
    nd={}
    for mask,val in dp.items():
        for col in range(6):
            if mask>>col&1 or not any(M[row][col]):continue
            key=mask|(1<<col);sign=-1 if (mask>>(col+1)).bit_count()%2 else 1
            prod=mul(val,M[row][col]);acc=nd.setdefault(key,[C.zero()]*3)
            for k in range(3):acc[k]+=sign*prod[k]
    dp=nd
    print('Sylvester row',row+1,'states',len(dp),'seconds',round(time.monotonic()-start,2),flush=True)
coeff=dp[63];assert not coeff[0] and not coeff[1]
expected=F[0].leading_coefficient()**3*G[0](0)*G[0](K((T-2)/r))*d['effective_second']
assert K(coeff[2])==expected
factor=str(coeff[2].factor());report={'status':'PASS_EXACT_RESULTANT_FIRST_NONZERO_COEFFICIENT',
 't_order':2,'t0_zero':True,'t1_zero':True,'root_formula_replay':True,
 'resultant_Hb_degree':118,'coefficient_factorization':factor,'elapsed_seconds':time.monotonic()-start,
 'scope':'Polynomial identity for all parameters. Nonvanishing requires old guards and T!=2, aq+1!=0, (q-1)T(aq+1)!=4. These three remaining walls are not excluded by this certificate.'}
save({'coefficients':coeff,'expected_root_formula':expected,'factorization':coeff[2].factor(),'series':d['series']},str(out/'pencil_resultant_initial_certificate.sobj'))
(out/'pencil_resultant_initial_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
