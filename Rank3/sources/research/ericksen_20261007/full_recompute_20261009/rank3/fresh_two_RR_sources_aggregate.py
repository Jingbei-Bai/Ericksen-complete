"""Assemble every fresh two-RR chart from independently computed source chunks."""
from sage.all import *
from pathlib import Path
import json,time,hashlib
out=Path(__file__).absolute().parent;st=time.monotonic();checks=[];matrices=[];reports=[]
def ck(n,v):
 if not v:raise AssertionError(n)
 checks.append(n)
for case in range(3):
 rows={};P=None
 for start in [0,6,12]:
  stem='case%d_rows%d_%d'%(case,start,start+6);p=out/'two_RR_chunks'/(stem+'.json');rec=json.loads(p.read_text());ck(stem+' passed',rec['status']=='PASS')
  for name,sha in rec['source_sha256'].items():ck(stem+' fresh input '+name,hashlib.sha256((out/name).read_bytes()).hexdigest()==sha)
  d=load(str(p.with_suffix('.sobj')));P=d['ring'];ck(stem+' exact source coverage',d['row_ids']==list(range(7))+list(range(7+start,7+start+6)))
  for j,index in enumerate(d['row_ids']):
   row=list(d['matrix'].row(j))
   if index in rows:ck(stem+' repeated initial source '+str(index),rows[index]==row)
   rows[index]=row
  for index in d['row_ids']:
   for col in range(5):ck(stem+' source verification record '+str((index,col)),'ALL twoRR source entries'+str((case,index,col)) in rec['checks'])
  reports.append(dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),source_object_sha256=hashlib.sha256(p.with_suffix('.sobj').read_bytes()).hexdigest(),seconds=rec['seconds']))
 ck('ALL25 original rows chart'+str(case),sorted(rows)==list(range(25)))
 matrices.append(matrix(P,[rows[i] for i in range(25)]))
save(dict(ring=P,matrices=matrices,source_chunks=reports),str(out/'fresh_two_RR_sources.sobj'),compress=False)
(out/'fresh_two_RR_sources.json').write_text(json.dumps(dict(status='PASS',count=len(checks),checks=checks,seconds=time.monotonic()-st,charts=3,entries=375,chunk_reports=reports,source='Each H50 is rebuilt from the fresh complete signed source DAG by exact homogeneous-product convolution. Full25x5 source matrices retain arbitrary symbolic a,b,H,Z,t,w,s,e. Source input hashes and repeated rows agree in all nine chunks; all375 coefficient comparisons appear in verified chunk records.',reproduction='Run fresh_two_RR_sources_chunked.py CASE START END for CASE=0,1,2 and (START,END)=(0,6),(6,12),(12,18), then fresh_two_RR_sources_aggregate.py.'),indent=2))
print('PASS all375 source entries;',len(checks),'assembly checks',time.monotonic()-st,flush=True)
