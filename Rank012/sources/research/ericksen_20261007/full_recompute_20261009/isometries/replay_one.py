"""Run a historical exact source in an isolated output tree, recording every load.

No saved PASS report is accepted as a substitute for executing source arithmetic.
Inputs produced by earlier fresh jobs override the immutable historical input.
"""
from pathlib import Path
import sys, json, hashlib, time, traceback, os, signal
from sage import all as sage

HERE = Path(__file__).absolute().parent
ORIGINAL = HERE.parents[1]
REPLAY = HERE/'replay'
rel = Path(sys.argv[1])
args = sys.argv[2:]
original = ORIGINAL/rel
if not original.is_file():
    original=HERE/'immutable_sources'/rel
target = REPLAY/rel
target.parent.mkdir(parents=True, exist_ok=True)
code = original.read_bytes().decode('utf8')
original_hash=hashlib.sha256(code.encode()).hexdigest()
target.write_text(code,encoding='utf8')
old_load = sage.load
loads = []
old_read_bytes = Path.read_bytes
old_read_text = Path.read_text
def fallback_path(p):
    p=p.absolute()
    if not p.exists() and p.is_relative_to(REPLAY):
        chosen=ORIGINAL/p.relative_to(REPLAY)
        if not chosen.is_file():chosen=HERE/'inputs'/p.relative_to(REPLAY)
        loads.append({'requested':str(p),'loaded':str(chosen),'provenance':'historical_exact_input',
                      'sha256':hashlib.sha256(old_read_bytes(chosen)).hexdigest()})
        return chosen
    return p
def tracked_bytes(p):return old_read_bytes(fallback_path(p))
def tracked_text(p,*a,**kw):return old_read_text(fallback_path(p),*a,**kw)
Path.read_bytes=tracked_bytes
Path.read_text=tracked_text
def tracked_load(filename, *a, **kw):
    p = Path(filename).absolute()
    chosen = p
    provenance = 'fresh_output'
    if not p.exists():
        chosen = ORIGINAL/p.relative_to(REPLAY)
        if not chosen.is_file():chosen=HERE/'inputs'/p.relative_to(REPLAY)
        provenance = 'historical_exact_input'
    loads.append({'requested':str(p),'loaded':str(chosen), 'provenance':provenance,
                  'sha256':hashlib.sha256(chosen.read_bytes()).hexdigest()})
    return old_load(str(chosen),*a,**kw)
sage.load = tracked_load
old_alarm=signal.alarm
signal.alarm=lambda seconds:old_alarm(3*seconds if seconds else 0)
sys.argv = [str(target)]+args
start = time.time()
status = 'PASS'
error = None
try:
    exec(compile(code,str(target),'exec'), {'__name__':'__main__','__file__':str(target)})
except BaseException as exc:
    status = 'FAIL'
    error = ''.join(traceback.format_exception(exc))
    print(error,flush=True)
finally:
    report={'status':status,'source':str(original),'source_sha256':original_hash,
            'executed_source_sha256':hashlib.sha256(code.encode()).hexdigest(),
            'source_text_transformations':[],
            'runtime_only_patch':'historical alarm budget multiplied by three; exact mathematics unchanged',
            'args':args,'seconds':time.time()-start,'loads':loads,'error':error}
    suffix='_'.join(args) if args else 'default'
    (target.parent/(target.stem+'_'+suffix+'_fresh_run.json')).write_text(json.dumps(report,indent=2),encoding='utf8')
if status != 'PASS':
    raise SystemExit(1)
