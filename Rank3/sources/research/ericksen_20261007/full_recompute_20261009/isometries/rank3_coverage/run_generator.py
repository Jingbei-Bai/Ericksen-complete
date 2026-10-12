"""Execute actual whole-source generators in an isolated, composable output tree."""
from pathlib import Path
import sys,json,hashlib,time,traceback,resource,os
from sage import all as sa
HERE=Path(__file__).absolute().parent
ROOT=HERE.parents[2]
ORIG=ROOT/'round10'
REPLAY=HERE/'replay'/'round10'
REPLAY.mkdir(parents=True,exist_ok=True)
name=sys.argv[1];args=sys.argv[2:];src=ORIG/(name+'.py');src=src if src.is_file() else HERE/'sources'/(name+'.py');target=REPLAY/src.name
code=src.read_bytes()
target.write_bytes(code)
loads=[];writes=[];old_load=sa.load;old_save=sa.save;old_text=Path.read_text;old_bytes=Path.read_bytes;old_cap=resource.setrlimit
old_resolve=Path.resolve
# The legacy per-process caps are execution choices, not mathematics.
resource.setrlimit=lambda *a,**kw:None
fresh_by_name={}
for q in (ROOT/'full_recompute_20261009'/'rank3'/'source_replays').glob('*/generated/round10/*.sobj'):
 record=q.parents[2]/'replay_record.json'
 if record.is_file() and json.loads(record.read_text()).get('status')=='PASS' and not (q.parents[2]/'execution.lock').exists():fresh_by_name[q.name]=q

def choose(p):
 p=Path(p).absolute()
 logical=ROOT/p.relative_to(HERE/'replay') if p.is_relative_to(HERE/'replay') else p
 views=json.loads(os.environ.get('ERICKSEN_INPUT_VIEWS','{}'))
 if logical.is_relative_to(ROOT.parent.parent):
  key=logical.relative_to(ROOT.parent.parent).as_posix()
  if key in views:
   mapped=ROOT.parent.parent/views[key]
   assert mapped.is_file(),str(mapped)
   return mapped,'fresh_declared_view'
 if p.exists():return p,'fresh_local' if p.is_relative_to(REPLAY) else 'direct_input'
 if p.is_relative_to(HERE/'replay'):
  if p.name in fresh_by_name:return fresh_by_name[p.name],'fresh_rank3_replay'
  return ROOT/p.relative_to(HERE/'replay'),'historical_coefficients_with_separate_source_bridge'
 return p,'direct_input'

def load(name,*a,**kw):
 p=Path(name);p=p if p.suffix else p.with_suffix('.sobj');chosen,prov=choose(p)
 loads.append(dict(path=str(chosen),sha256=hashlib.sha256(old_bytes(chosen)).hexdigest(),provenance=prov))
 return old_load(str(chosen),*a,**kw)

def save(value,name,*a,**kw):
 pp=Path(name).absolute()
 assert pp.is_relative_to(HERE),str(pp)
 writes.append(str(pp));return old_save(value,str(pp),*a,**kw)

def text(p,*a,**kw):
 chosen,prov=choose(p)
 if chosen!=p:loads.append(dict(path=str(chosen),sha256=hashlib.sha256(old_bytes(chosen)).hexdigest(),provenance=prov))
 return old_text(chosen,*a,**kw)
sa.load=load;sa.save=save;Path.read_text=text
Path.resolve=lambda self,*args,**kwargs:self.absolute()
sys.argv=[str(target)]+args
ns={'__name__':'__main__','__file__':str(target)}
start=time.monotonic();status='PASS';error=None;extras=[]
try:
 try:
  exec(compile(code,str(target),'exec'),ns)
 except SystemExit as exc:
  if exc.code not in (None,0):raise
 if name in ['independent_subcritical_single_RR_complete_E_relative_dual','subcritical_single_full_conic_G2_relative_dual']:
  assert ns['ans']=={},'Fresh whole nonreduced source remainder is not empty'
  extras.append('fresh nonreduced conic derivative remainder exactly zero')
except BaseException as exc:
 status='FAIL';error=traceback.format_exc();print(error,flush=True)
finally:
 sa.load=old_load;sa.save=old_save;Path.read_text=old_text;Path.resolve=old_resolve;resource.setrlimit=old_cap
 rec=dict(status=status,source=str(src),source_sha256=hashlib.sha256(old_bytes(src)).hexdigest(),executed_sha256=hashlib.sha256(code).hexdigest(),source_text_transformations=[],args=args,seconds=time.monotonic()-start,loads=loads,writes=writes,error=error,additional_checks=extras)
 (HERE/(name+('_'+'_'.join(args) if args else '')+'_run.json')).write_text(json.dumps(rec,indent=2))
 print('FRESH WHOLE SOURCE',name,args,status,rec['seconds'],flush=True)
if status!='PASS':sys.exit(1)
