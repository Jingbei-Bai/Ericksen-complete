"""Execute a whole source with the exact ray comparison metadata adapter."""
from pathlib import Path
import hashlib
import json
import sys
import time
import traceback
from sage import all as sa
from reproduction.ray_fixture import restore_ray_fixture

ROOT = Path(__file__).absolute().parent.parent
member = Path(sys.argv[1])
if member.is_absolute():
    source = member.absolute()
    assert source.is_relative_to(ROOT), source
else:
    assert '..' not in member.parts
    source = ROOT / member
code = source.read_bytes()
digest = hashlib.sha256(code).hexdigest()
receipt = ROOT / 'reproduction/reports' / (source.stem + '_execution.json')
receipt.parent.mkdir(parents=True, exist_ok=True)
start = time.monotonic()
load0 = sa.load
reads = []


def load(name, *args, **kwargs):
    path = Path(name)
    if not path.suffix:
        path = path.with_suffix('.sobj')
    value = load0(str(path), *args, **kwargs)
    reads.append({'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                  'fixture_kind': value.get('fixture_kind') if isinstance(value, dict) else None})
    return restore_ray_fixture(value)


sa.load = load
status, error = 'PASS', None
sys.argv = [str(source)] + sys.argv[2:]
try:
    exec(compile(code, str(source), 'exec'), {'__name__': '__main__', '__file__': str(source)})
except BaseException as exc:
    status, error = 'FAIL', repr(exc)
    traceback.print_exc()
finally:
    sa.load = load0
    assert hashlib.sha256(source.read_bytes()).hexdigest() == digest
    receipt.write_text(json.dumps({'status': status, 'source': str(source),
                                  'executed_source_sha256': digest, 'source_text_transformations': [],
                                  'imports': reads, 'error': error, 'seconds': time.monotonic() - start},
                                 indent=2) + '\n')
    print('RAY METADATA FULL SOURCE EXECUTION', status, source.stem, digest, flush=True)
if status != 'PASS':
    raise SystemExit(1)
