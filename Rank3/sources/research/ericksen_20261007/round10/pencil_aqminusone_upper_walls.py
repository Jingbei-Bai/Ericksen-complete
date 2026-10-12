"""Higher bounded determinant coefficients, then restriction to the two walls."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'pencil_aqminusone_newton.sobj'))
C=PolynomialRing(QQ,names=('q','T'));q,T=C.gens();N=9
vectors=[]
for p,ht in zip(d['cubics'],(18,21)):
    vectors.append([[C({tuple(mon[:2]):val for mon,val in co.dict().items() if mon[2]==ht-k}) for k in range(N)] for co in reversed(p.list())])
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
    dp=nd;print('row',row,'sec',round(time.monotonic()-start,2),flush=True)
coeff=dp[63]
save({'coefficients':coeff,'vectors':vectors},str(out/'pencil_aqminusone_upper_walls.sobj'))
print('SAVED',round(time.monotonic()-start,2),flush=True)
R=PolynomialRing(QQ,'q');qq=R.gen();K=R.fraction_field();S=PolynomialRing(K,'T');tt=S.gen()
phi=C.hom([S(qq),tt],S);NN=(qq*qq+6*qq+1)*tt-4*(qq-1)**2;JJ=(tt-1)*NN-8*qq*tt
records=[];remainders={}
for name,wall in [('N',NN),('J',JJ)]:
    remainders[name]=[]
    for k in range(4,N):
        rem=phi(coeff[k])%wall;remainders[name].append(rem)
        rec={'wall':name,'t_order':k,'zero':not bool(rem)}
        if rem:
            rec['factorization']=str(rem.factor());records.append(rec);print(rec,flush=True);break
        records.append(rec);print(rec,flush=True)
save({'coefficients':coeff,'vectors':vectors,'remainders':remainders,'records':records,'walls':(NN,JJ)},str(out/'pencil_aqminusone_upper_walls.sobj'))
(out/'pencil_aqminusone_upper_walls.json').write_text(json.dumps({'status':'EXACT_HIGHER_COEFFICIENT_WALL_REMAINDERS','records':records,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
