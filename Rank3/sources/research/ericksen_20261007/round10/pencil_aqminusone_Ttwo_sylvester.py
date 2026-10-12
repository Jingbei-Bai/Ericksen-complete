"""Independent finite determinant check of the T=2 root-lifting defect."""
from pathlib import Path
import json,time
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'pencil_aqminusone_Ttwo_defect.sobj'))
F,G=d['series'];C=F[0].base_ring().ring();q=C.gen();N=3
vectors=[[[C(row[k][j]) for k in range(N)] for j in reversed(range(4))] for row in (F,G)]
zero=[C.zero()]*N;M=[]
for vv in vectors:
    for shift in range(3):M.append([zero]*shift+vv+[zero]*(2-shift))
def mul(x,y):return [sum((x[j]*y[k-j] for j in range(k+1)),C.zero()) for k in range(N)]
dp={0:[C.one()]+[C.zero()]*(N-1)}
for row in range(6):
    nd={}
    for mask,val in dp.items():
        for col in range(6):
            if mask>>col&1 or not any(M[row][col]):continue
            key=mask|(1<<col);sign=-1 if (mask>>(col+1)).bit_count()%2 else 1
            prod=mul(val,M[row][col]);acc=nd.setdefault(key,[C.zero()]*N)
            for k in range(N):acc[k]+=sign*prod[k]
    dp=nd
coeff=dp[63];assert not coeff[0] and not coeff[1] and coeff[2]
report={'status':'PASS_EXACT_TTWO_INTERSECTION','Hb_resultant_power':100,'coefficient':str(coeff[2].factor()),'lower_t_coefficients_zero':True,'elapsed_seconds':time.monotonic()-start}
save({'coefficients':coeff,'vectors':vectors},str(out/'pencil_aqminusone_Ttwo_sylvester.sobj'))
(out/'pencil_aqminusone_Ttwo_sylvester.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
