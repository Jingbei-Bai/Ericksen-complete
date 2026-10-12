"""The aq=-1,T=2 intersection, exact single-parameter Newton data."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'pencil_aqminusone_newton.sobj'))
C=PolynomialRing(QQ,names=('q','Hb'));q,Hb=C.gens();S=PolynomialRing(C,'U');U=S.gen();polys=[];records=[]
for i,p in enumerate(d['cubics']):
    cc=[];rr=[]
    for j,co in enumerate(p.list()):
        acc={}
        for mon,val in co.dict().items():
            key=(mon[0],mon[2]);acc[key]=acc.get(key,0)+val*2**mon[1]
        f=C(acc);cc.append(f)
        lo=min(m[1] for m in f.dict());hi=max(m[1] for m in f.dict())
        def at(n):return C({(m[0],0):v for m,v in f.dict().items() if m[1]==n})
        rec={'power_U':j,'degree_Hb':int(hi),'valuation_Hb':int(lo),'leading':str(at(hi).factor()),'lowest':str(at(lo).factor())};rr.append(rec);print(i,rec,flush=True)
    polys.append(S(cc));records.append(rr)
save({'base_ring':C,'cubics':polys,'records':records},str(out/'pencil_aqminusone_Ttwo.sobj'))
(out/'pencil_aqminusone_Ttwo.json').write_text(json.dumps({'records':records,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
