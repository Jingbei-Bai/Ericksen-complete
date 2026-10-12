"""Exact specialization aq=-1; both remaining parameters q,T stay symbolic."""
from pathlib import Path
import json,time
from sage.all import load,save,QQ,PolynomialRing
out=Path(__file__).resolve().parent;start=time.monotonic()
d=load(str(out/'pencil_vzero_sone_cubics.sobj'))
C=PolynomialRing(QQ,names=('q','T','Hb'));q,T,Hb=C.gens()
S=PolynomialRing(C,names=('U',));U=S.gen();polys=[];recs=[]
for idx,p in enumerate(d['cubics']):
    # A common q power is harmless on this wall: a=-1/q>0 implies q<0.
    e=max(int(mon[0]-mon[1]) for co in p.list() for mon in co.dict())
    cc=[];record=[]
    for j,co in enumerate(p.list()):
        acc={}
        for mon,val in co.dict().items():
            key=(int(mon[1]-mon[0]+e),int(mon[2]),int(mon[3]))
            acc[key]=acc.get(key,0)+(-1)**int(mon[0])*val
        f=C(acc);cc.append(f)
        minH=min(mon[2] for mon in f.dict());maxH=max(mon[2] for mon in f.dict())
        def at(n):return C({(mon[0],mon[1],0):v for mon,v in f.dict().items() if mon[2]==n})
        rec={'U_power':j,'Hb_valuation':int(minH),'Hb_degree':int(maxH),'lowest':str(at(minH).factor()),'highest':str(at(maxH).factor()),'terms':f.number_of_terms()}
        record.append(rec);print('row',idx,rec,flush=True)
    f=S(cc);polys.append(f);recs.append({'source':idx,'q_clearing_power':e,'coefficients':record})
save({'base_ring':C,'slice_ring':S,'cubics':polys,'records':recs,'specialization':'a=-1/q; specialized cubic multiplied by q^e, e in records'},str(out/'pencil_aqminusone_newton.sobj'))
(out/'pencil_aqminusone_newton.json').write_text(json.dumps({'status':'EXACT_WALL_SPECIALIZATION','records':recs,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
print('DONE',round(time.monotonic()-start,2),flush=True)
