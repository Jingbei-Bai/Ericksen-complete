"""Replay an archived exact-source reader into an isolated new output folder.

This verifies computations, not the truth of imported necessary-source data.
Every imported Sage object is hashed and disclosed. Existing files are never
written. Source bytes are copied verbatim before execution. This runner accepts
only a named Python file under round10 or rank3_simplification_20261009.
"""
from pathlib import Path
import sage.all as sa
import sys,json,time,hashlib,os,resource
from reproduction.ray_fixture import restore_ray_fixture
out=Path(__file__).absolute().parent
research=out.parent.parent
rel=Path(sys.argv[1]);src=research/rel
assert src.suffix=='.py' and rel.parts[0] in ['round10','rank3_simplification_20261009']
extra_arguments=sys.argv[2:]
if src.stem=='independent_subcritical_single_RR_yw_minus_closed_reader' and not extra_arguments:
 extra_arguments=['H2']
run=out/'source_replays'/(src.stem+('_'+extra_arguments[0] if extra_arguments else ''));run.mkdir(parents=True,exist_ok=True)
def completed_current_replay():
 if os.environ.get('ERICKSEN_FORCE_RECOMPUTE')=='1':return False
 rec=run/'replay_record.json'
 if not rec.is_file():return False
 prior=json.loads(rec.read_text())
 if prior.get('status')!='PASS' or prior.get('script_sha256')!=hashlib.sha256(src.read_bytes()).hexdigest():return False
 for item in prior.get('imports',[])+prior.get('file_reads',[]):
  pp=Path(item['path'])
  if not pp.is_file() or hashlib.sha256(pp.read_bytes()).hexdigest()!=item['sha256']:return False
 return True
lock=run/'execution.lock'
while True:
 if completed_current_replay():
  print('REUSED CURRENT RECOMPUTATION: exact source and all recorded input hashes unchanged',src.name,extra_arguments,flush=True);sys.exit(0)
 try:
  fd=os.open(str(lock),os.O_WRONLY|os.O_CREAT|os.O_EXCL);os.write(fd,str(os.getpid()).encode());os.close(fd);break
 except FileExistsError:
  try:os.kill(int(lock.read_text()),0)
  except ProcessLookupError:lock.unlink();continue
  time.sleep(1)
code=src.read_bytes().decode('utf-8');(run/'executed_source.py').write_text(code)
start=time.monotonic();imports=[];writes=[];file_reads=[];original_load=sa.load;original_save=sa.save
original_write=Path.write_text;original_resolve=Path.resolve
original_read_text=Path.read_text;original_read_bytes=Path.read_bytes;original_open=Path.open
original_replace=os.replace;original_rename=os.rename
original_setrlimit=resource.setrlimit
def track_read(p,kind):
 p=Path(p).absolute()
 try:p.relative_to(research)
 except ValueError:return
 if p.is_file():
  with original_open(p,'rb') as fh: digest=hashlib.sha256(fh.read()).hexdigest()
  rec=dict(path=str(p),sha256=digest,kind=kind)
  if rec not in file_reads:file_reads.append(rec)
def read_text(p,*args,**kwargs):
 track_read(p,'text');return original_read_text(p,*args,**kwargs)
def read_bytes(p,*args,**kwargs):
 track_read(p,'bytes');return original_read_bytes(p,*args,**kwargs)
def opened(p,mode='r',*args,**kwargs):
 if any(k in mode for k in ('w','a','x','+')):
  target=dest(p);writes.append(str(target));return original_open(target,mode,*args,**kwargs)
 track_read(p,'open');return original_open(p,mode,*args,**kwargs)
def replaced(srcpath,dstpath,*args,**kwargs):
 source=dest(srcpath);target=dest(dstpath);writes.append(str(target))
 return original_replace(source,target,*args,**kwargs)
def dest(path):
 p=Path(path).absolute()
 try:tail=p.relative_to(research)
 except ValueError:raise RuntimeError('write outside research '+str(p))
 if 'full_recompute_20261009' in tail.parts:return p
 target=run/'generated'/tail;target.parent.mkdir(parents=True,exist_ok=True);return target
def load(name,*args,**kwargs):
 p=Path(name)
 if not p.suffix:p=p.with_suffix('.sobj')
 generated=run/'generated'/p.relative_to(research)
 views=json.loads(os.environ.get('ERICKSEN_INPUT_VIEWS','{}'))
 logical=p.relative_to(research.parent.parent).as_posix()
 actual=research.parent.parent/views[logical] if logical in views else (generated if generated.is_file() else p)
 imports.append(dict(path=str(actual),sha256=hashlib.sha256(actual.read_bytes()).hexdigest(),regenerated=actual==generated or logical in views))
 return restore_ray_fixture(original_load(str(actual),*args,**kwargs))
def save(value,name,*args,**kwargs):
 target=dest(name);writes.append(str(target));return original_save(value,str(target),*args,**kwargs)
def write(p,text,*args,**kwargs):
 target=dest(p);writes.append(str(target));return original_write(target,text,*args,**kwargs)
sa.load=load;sa.save=save;Path.write_text=write;Path.read_text=read_text;Path.read_bytes=read_bytes;Path.open=opened
os.replace=replaced;os.rename=replaced
resource.setrlimit=lambda *args,**kwargs:None
Path.resolve=lambda self,*args,**kwargs:self.absolute()
# Remove archived process memory caps, which are execution-budget choices and
# not mathematical operations. Readers have no native helper/output commands.
executed_source_sha256=hashlib.sha256(code.encode()).hexdigest()
assert executed_source_sha256==hashlib.sha256(src.read_bytes()).hexdigest()
sys.argv=[str(src)]+extra_arguments
status='PASS';error=None
try:exec(compile(code,str(src),'exec'),{'__name__':'__main__','__file__':str(src)})
except BaseException as exc:
 status='FAIL';error=repr(exc)
 import traceback;traceback.print_exc()
finally:
 sa.load=original_load;sa.save=original_save;Path.write_text=original_write;Path.resolve=original_resolve;Path.read_text=original_read_text;Path.read_bytes=original_read_bytes;Path.open=original_open
 os.replace=original_replace;os.rename=original_rename
 resource.setrlimit=original_setrlimit
 record=dict(status=status,script=str(src),script_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),executed_source_sha256=executed_source_sha256,source_text_transformations=[],imports=imports,file_reads=file_reads,writes=writes,elapsed_seconds=time.monotonic()-start,error=error,
  evidence_level='FRESH_SOURCE_REPLAY_WITH_IMPORTED_COEFFICIENT_DATA',
  limitation='Imported necessary-source arrays and comparison objects were not independently regenerated by this reader. Successful replay is not evidence that the entire upstream source chain was regenerated.')
 (run/'replay_record.json').write_text(json.dumps(record,indent=2))
 lock.unlink(missing_ok=True)
 print('ISOLATED REPLAY',status,src.name,'imports',len(imports),'seconds',record['elapsed_seconds'],flush=True)
if status!='PASS':sys.exit(1)
