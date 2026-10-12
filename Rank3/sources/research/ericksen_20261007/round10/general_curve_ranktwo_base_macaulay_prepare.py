"""Exact source matrix plus bounded modular support discovery for a QQ(H) lift."""
from sage.all import *
from pathlib import Path
from cysignals.alarm import alarm,cancel_alarm
from cysignals.signals import AlarmInterrupt
import json,time
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'general_curve_canonical_ranktwo_linear.sobj'));C=d['ring'];Hc,Zc,Sc=C.gens()
Hp=PolynomialRing(QQ,'H');H=Hp.gen();F=Hp.fraction_field();R=PolynomialRing(F,names=('Z','S'));Z,S=R.gens()
def cv(f):
    terms={}
    for e,c in C(f).dict().items():
        ee=(int(e[1]),int(e[2]));terms[ee]=terms.get(ee,Hp.zero())+c*H**int(e[0])
    return R(terms)
sources=[cv(f) for f in d['parameter_sources']];target=(H-Z)**2*(Z-1)**2*S**2*(S-1)**2;ff=GF(32003);h0=ff(2)
def ev(c):return ff(c.numerator()(2))/ff(c.denominator()(2))
records=[];success=False
alarm(25)
try:
    for deg in range(22,26):
        mons=[(i,n-i) for n in range(deg+1) for i in range(n+1)];idx={e:i for i,e in enumerate(mons)};prov=[];column_dicts=[]
        for j,f in enumerate(sources):
            for n in range(deg-int(f.total_degree())+1):
                for k in range(n+1):
                    exp=(k,n-k);prov.append((j,exp));column_dicts.append({(int(e[0])+exp[0],int(e[1])+exp[1]):c for e,c in f.dict().items()})
        mm=matrix(ff,len(mons),len(prov),sparse=False)
        for j,col in enumerate(column_dicts):
            for e,c in col.items():mm[idx[e],j]=ev(c)
        bv=vector(ff,[ev(target.monomial_coefficient(Z**e[0]*S**e[1])) for e in mons]);cc=mm.pivots();rr=mm.matrix_from_columns(cc).transpose().pivots();x=mm.matrix_from_rows_and_columns(rr,cc).solve_right(vector(ff,[bv[i] for i in rr]));ok=mm.matrix_from_columns(cc)*x==bv
        rec={'degree':deg,'rows':len(mons),'cols':len(prov),'rank':len(cc),'target_member_modular':bool(ok),'seconds':time.monotonic()-start};records.append(rec);print(rec,flush=True)
        if ok:
            keep=[i for i in range(len(cc)) if x[i]];cols=[cc[i] for i in keep];rows=mm.matrix_from_columns(cols).transpose().pivots();assert len(rows)==len(cols);success=True;break
except (AlarmInterrupt,KeyboardInterrupt):pass
finally:cancel_alarm()
payload={'ring':R,'sources':sources,'target':target,'records':records,'scope':'Only the source matrix is exact. Selected pivots and modular membership are discovery data, not an identity.'}
if success:
    A=matrix(F,len(mons),len(cols))
    for j,k in enumerate(cols):
        for e,c in column_dicts[k].items():A[idx[e],j]=c
    # Every retained source column and the complete target are reconstructed.
    assert all(sum(A[i,j]*Z**mons[i][0]*S**mons[i][1] for i in range(len(mons)))==sources[prov[k][0]]*Z**prov[k][1][0]*S**prov[k][1][1] for j,k in enumerate(cols))
    b=vector(F,[target.monomial_coefficient(Z**e[0]*S**e[1]) for e in mons]);assert sum(b[i]*Z**e[0]*S**e[1] for i,e in enumerate(mons))==target
    payload.update(matrix=A,target_vector=b,monomials=mons,provenance=[prov[k] for k in cols],selected_rows=rows,all_source_roundtrips=True,target_roundtrip=True)
tmp=out/'general_curve_ranktwo_base_macaulay_input_next.sobj';save(payload,str(tmp));tmp.replace(out/'general_curve_ranktwo_base_macaulay_input.sobj')
rep={'modular_support_found':success,'records':records,'selected_size':len(cols) if success else None,'seconds':time.monotonic()-start,'scope':payload['scope']};(out/'general_curve_ranktwo_base_macaulay_prepare.json').write_text(json.dumps(rep,indent=2)+'\n');print(rep,flush=True)
