"""Factored coefficient norms for the two residual Newton walls."""
from pathlib import Path
import json,time
from sage.all import load,save
out=Path(__file__).resolve().parent;start=time.monotonic();d=load(str(out/'pencil_aqminusone_upper_walls.sobj'))
records=[];data={}
for name,wall in zip(('N','J'),d['walls']):
 rem=next(p for p in d['remainders'][name] if p);norm=wall.resultant(rem)
 facnum=norm.numerator().factor();facden=norm.denominator().factor()
 rec={'wall':name,'norm_numerator':str(facnum),'norm_denominator':str(facden)};print(rec,flush=True);records.append(rec)
 data[name]={'remainder':rem,'norm':norm,'norm_factorization':(facnum,facden)}
save(data,str(out/'pencil_aqminusone_wall_factors.sobj'))
(out/'pencil_aqminusone_wall_factors.json').write_text(json.dumps({'records':records,'elapsed_seconds':time.monotonic()-start},indent=2)+'\n')
