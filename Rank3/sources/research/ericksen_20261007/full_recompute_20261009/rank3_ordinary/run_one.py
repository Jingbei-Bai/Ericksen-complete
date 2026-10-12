"""Isolated algorithm replay. Rebuilt source overrides cached witnesses."""
from pathlib import Path
import sage.all as sa
import json,time,sys,hashlib,traceback
import resource,os
from reproduction.ray_fixture import restore_ray_fixture
# Historical small address-space limits are execution settings, not mathematics.
resource.setrlimit=lambda *args,**kwargs: None
out=Path(__file__).absolute().parent;root=out.parent.parent;mirror=out/'rebuild'
rel=Path(sys.argv[1]);extra_args=sys.argv[2:];src=root/rel;target=mirror/rel;target.parent.mkdir(parents=True,exist_ok=True)
code=src.read_bytes().decode('utf-8');target.write_text(code);start=time.monotonic();imports=[]
load0=sa.load;read0=Path.read_bytes;readtext0=Path.read_text;open0=Path.open
def locate(p):
 p=Path(p).absolute()
 if p.exists():return p
 try:q=root/p.relative_to(mirror)
 except ValueError:return p
 return q if q.exists() else p
def record(p,kind):
 if p.is_file():imports.append(dict(path=str(p),regenerated=p.is_relative_to(mirror),sha256=hashlib.sha256(read0(p)).hexdigest(),access=kind))
def rd(p,*a,**k):
 actual=locate(p);record(actual,'read_bytes');return read0(actual,*a,**k)
def rt(p,*a,**k):
 actual=locate(p);record(actual,'read_text');return readtext0(actual,*a,**k)
def op(p,mode='r',*a,**k):
 if 'r' in mode and not any(v in mode for v in 'wax+'):
  p=locate(p)
  # Avoid recursion through read_bytes when hashing an open call.
  if p.is_file():
   with open0(p,'rb') as f:sha=hashlib.sha256(f.read()).hexdigest()
   imports.append(dict(path=str(p),regenerated=p.is_relative_to(mirror),sha256=sha,access='open'))
 return open0(p,mode,*a,**k)
def ld(name,*a,**k):
 p=Path(name)
 if not p.suffix:p=p.with_suffix('.sobj')
 actual=locate(p);imports.append(dict(path=str(actual),regenerated=actual.is_relative_to(mirror),sha256=hashlib.sha256(read0(actual)).hexdigest()))
 return restore_ray_fixture(load0(str(actual),*a,**k))
sa.load=ld;Path.read_bytes=rd;Path.read_text=rt;Path.open=op;Path.resolve=lambda self,*a,**k:self.absolute()
status='PASS';error=None
executed_source_sha256=hashlib.sha256(code.encode()).hexdigest()
assert executed_source_sha256==hashlib.sha256(src.read_bytes()).hexdigest()
sys.argv=[str(target)]+sys.argv[2:]
target.write_text(code)
executed_source_sha256=hashlib.sha256(code.encode()).hexdigest()
try:exec(compile(code,str(target),'exec'),dict(__name__='__main__',__file__=str(target)))
except BaseException as exc:status='FAIL';error=repr(exc);traceback.print_exc()
finally:
 sa.load=load0;Path.read_bytes=read0;Path.read_text=readtext0;Path.open=open0
 report=dict(status=status,script=str(rel),arguments=extra_args,imports=imports,error=error,seconds=time.monotonic()-start,executed_source_sha256=executed_source_sha256,source_text_transformations=[])
 (out/'job_reports').mkdir(exist_ok=True);(out/'job_reports'/(src.stem+('_'+extra_args[0] if extra_args else '')+'.json')).write_text(json.dumps(report,indent=2));print('JOB',status,rel,report['seconds'],flush=True)
if status!='PASS':sys.exit(1)
