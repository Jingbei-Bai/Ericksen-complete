"""Exact all-parameter highest-u coefficient needed for localization at s."""
from pathlib import Path
import time,json
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_base_content_and_klein.sobj'));R=d['fiber_ring'];C=d['base_ring'];a,q,T,Hb=C.gens()
p=d['primitive_rows'][0];uu,vv,ss=R.gens();top={tuple(mon):co for mon,co in p.dict().items() if mon[0]==p.degree(uu)}
assert set(top)=={(6,0,1)},top
c6=C(top[(6,0,1)])
facts=list(c6.factor());records=[];witnesses=[]
for f,e in facts:
    hd=f.degree(Hb);lead=C({tuple(0 if j==3 else ee for j,ee in enumerate(mm)):cc for mm,cc in f.dict().items() if mm[3]==hd})
    rec={'factor':str(f),'power':int(e),'Hb_degree':int(hd),'Hb_leading_coefficient_factorization':str(lead.factor())}
    records.append(rec);witnesses.append((f,e,lead))
save({'coefficient':c6,'factors':facts,'Hb_leading_witnesses':witnesses,'source_identity':p.coefficient({uu:6})==R(c6)*ss},str(out/'pencil_u_primitive_coefficient.sobj'))
report={'status':'EXACT_U_LEADING_COEFFICIENT','u_degree':6,'only_top_monomial':[6,0,1],'factors':records,'elapsed_seconds':time.monotonic()-start,
        'scope':'Nonzero c6 makes H primitive in u over the coefficient ring localized at s. Uniform parameter nonvanishing still requires examining all factor coefficients.'}
(out/'pencil_u_primitive_coefficient.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2),flush=True)
