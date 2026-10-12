"""Select a source lift for B*S*(S-1) using only the three determinant coefficients."""
from sage.all import *
from pathlib import Path
from cysignals.alarm import alarm,cancel_alarm
from cysignals.signals import AlarmInterrupt
import json,time
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'general_curve_zero_det_coefficients_input.sobj'));C=d['ring'];Hc,Zc,Sc=C.gens();Hp=PolynomialRing(QQ,'H');H=Hp.gen();F=Hp.fraction_field();R=PolynomialRing(F,names=('Z','S'));Z,S=R.gens()
def cv(f):
    terms={}
    for e,c in C(f).dict().items():ee=(int(e[1]),int(e[2]));terms[ee]=terms.get(ee,Hp.zero())+c*H**int(e[0])
    return R(terms)
src=[cv(f) for f in d['sources']];B=H*Z*(H**2-1)*S**2+2*(Z**2-H**2)*S-H*Z*(Z**2-1);target=B*S*(S-1);ff=GF(32003)
def ev(c,h=3):return ff(c.numerator()(h))/ff(c.denominator()(h))
records=[];payload={'ring':R,'sources':src,'target':target,'source_file':'general_curve_zero_det_coefficients_input.sobj','scope':'Modular selection only; every selected column and target is reconstructed over QQ(H), and a later exact lift is required.'};print('SOURCE_DEGREES',[int(f.total_degree()) for f in src],flush=True);alarm(25)
try:
    for D in range(max(int(f.total_degree()) for f in src),35):
        mons=[(k,n-k) for n in range(D+1) for k in range(n+1)];idx={e:i for i,e in enumerate(mons)};cols=[];prov=[]
        for j,f in enumerate(src):
            for n in range(D-int(f.total_degree())+1):
                for k in range(n+1):ex=(k,n-k);cols.append({(int(e[0])+ex[0],int(e[1])+ex[1]):c for e,c in f.dict().items()});prov.append((j,ex))
        M=matrix(ff,len(mons),len(cols))
        for j,col in enumerate(cols):
            for e,c in col.items():M[idx[e],j]=ev(c)
        b=vector(ff,[ev(target.monomial_coefficient(Z**e[0]*S**e[1])) for e in mons]);cc=M.pivots();rr=M.matrix_from_columns(cc).transpose().pivots();x=M.matrix_from_rows_and_columns(rr,cc).solve_right(vector(ff,[b[i] for i in rr]));ok=M.matrix_from_columns(cc)*x==b;rec={'degree':D,'shape':[len(mons),len(cols)],'rank':len(cc),'modular_member':bool(ok),'seconds':time.monotonic()-start};records.append(rec);print(rec,flush=True)
        if ok:
            cc=[cc[j] for j in range(len(cc)) if x[j]];rr=M.matrix_from_columns(cc).transpose().pivots();A=matrix(F,len(mons),len(cc));bv=vector(F,[target.monomial_coefficient(Z**e[0]*S**e[1]) for e in mons])
            for j,k in enumerate(cc):
                for e,c in cols[k].items():A[idx[e],j]=c
            assert all(sum(A[i,j]*Z**e[0]*S**e[1] for i,e in enumerate(mons))==src[prov[k][0]]*Z**prov[k][1][0]*S**prov[k][1][1] for j,k in enumerate(cc));assert sum(bv[i]*Z**e[0]*S**e[1] for i,e in enumerate(mons))==target
            stability=[]
            for h in [2,5,7]:
                mm=matrix(ff,[[ev(c,h) for c in row] for row in A.rows()]);bb=vector(ff,[ev(c,h) for c in bv]);ss=mm.matrix_from_rows(rr);rk=ss.rank();it={'H':h,'rank':int(rk)}
                if rk==A.ncols():xx=ss.solve_right(vector(ff,[bb[i] for i in rr]));it['all_rows_match']=mm*xx==bb
                stability.append(it)
            payload.update(matrix=A,target_vector=bv,selected_rows=rr,monomials=mons,provenance=[prov[k] for k in cc],source_and_target_roundtrips=True,stability=stability,status='EXACT_INPUT_PREPARED_MODULAR_SUPPORT_ONLY');break
    else:payload['status']='NO_SUPPORT_FOUND'
except (AlarmInterrupt,KeyboardInterrupt):payload['status']='BOUNDED_UNKNOWN'
finally:cancel_alarm()
payload['records']=records;save(payload,str(out/'general_curve_zero_det_B_macaulay_input.sobj'));rep={'status':payload['status'],'records':records,'selected_size':payload['matrix'].ncols() if 'matrix' in payload else None,'stability':payload.get('stability'),'seconds':time.monotonic()-start};(out/'general_curve_zero_det_B_macaulay_prepare.json').write_text(json.dumps(rep,indent=2,default=int)+'\n');print(rep,flush=True)
