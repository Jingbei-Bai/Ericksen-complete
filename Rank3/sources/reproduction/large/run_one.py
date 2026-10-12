"""Execute a large calculation source in its newly generated input directory."""
from pathlib import Path
import hashlib, json, os, sys, tempfile, traceback
import sage.all as sage

out=Path(__file__).resolve().parent
root=out.parent.parent
relative=Path(sys.argv[1]);extra=sys.argv[2:]
assert relative.parts[0]=='round10' and relative.suffix=='.py'
original=root/'research/ericksen_20261007'/relative
target=out/'generated'/relative;target.parent.mkdir(parents=True,exist_ok=True)
def publish_source(source_path, destination):
    content=source_path.read_bytes()
    if destination.is_file() and destination.read_bytes()==content:
        return content
    with tempfile.NamedTemporaryFile(dir=destination.parent, prefix='.source-', delete=False) as handle:
        temporary=Path(handle.name)
        handle.write(content)
    os.replace(temporary,destination)
    assert destination.read_bytes()==content
    return content

source=publish_source(original,target)
if original.name=='general_curve_zero_det_B_interpolation_prepare.py':
    sibling=original.with_name('general_curve_ranktwo_integer_interpolation.py')
    publish_source(sibling,target.with_name(sibling.name))
digest=hashlib.sha256(source).hexdigest()
sys.argv=[str(target)]+extra
status='PASS';error=None
inputs=[];load_original=sage.load
def load_bound(path,*arguments,**keywords):
    resolved=Path(str(path)).resolve()
    if not resolved.is_file() and Path(str(resolved)+'.sobj').is_file():
        resolved=Path(str(resolved)+'.sobj')
    assert resolved.is_file(),str(resolved)
    before=hashlib.sha256(resolved.read_bytes()).hexdigest()
    value=load_original(str(resolved),*arguments,**keywords)
    assert hashlib.sha256(resolved.read_bytes()).hexdigest()==before
    inputs.append({'path':str(resolved),'sha256':before,
                   'generated':resolved.is_relative_to(out/'generated')})
    return value
sage.load=load_bound
try:
    exec(compile(source,str(target),'exec'),{'__name__':'__main__','__file__':str(target)})
except BaseException as exception:
    status='FAIL';error=repr(exception);traceback.print_exc()
finally:
    sage.load=load_original
report={'status':status,'executed_source_sha256':digest,'source_text_transformations':[],
    'original_source':str(original),'executed_source':str(target),'imports':inputs,'error':error}
records=out/'job_reports';records.mkdir(exist_ok=True)
(records/(original.stem+('_'+extra[0] if extra else '')+'.json')).write_text(json.dumps(report,indent=2)+'\n')
if status!='PASS': raise SystemExit(1)
