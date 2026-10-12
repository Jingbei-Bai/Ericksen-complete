"""Accept current-run results only while their algebraic inputs are unchanged."""
from pathlib import Path
import hashlib,json
OUT=Path(__file__).absolute().parent
RESEARCH=OUT.parent.parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def local(value):
 if value.startswith('/mnt/c/') and OUT.drive:value='C:/'+value[7:]
 return Path(value)
def valid_report(path):
 if not path.exists():return False
 d=json.loads(path.read_text())
 if d.get('status')!='PASS':return False
 stem=Path(d['script']).stem
 executed=OUT/'rebuild'/d['script']
 if d.get('executed_source_sha256') and (not executed.is_file() or sha(executed)!=d['executed_source_sha256']):return False
 for item in d['imports']:
  if item['path']=='/proc/self/statm':continue
  p=local(item['path'])
  # A producer may read its own atomic output as a resumable checkpoint.
  if p.name.startswith(stem+'.') and p.suffix=='.sobj':continue
  if p.is_relative_to(RESEARCH) and not p.is_relative_to(OUT):
   preferred=OUT/'rebuild'/p.relative_to(RESEARCH)
   if preferred.is_file() and sha(preferred)!=item['sha256']:
    # Exact comparison may justify a newer serialization of the same
    # input; corrected objects never qualify for this bridge.
    comparison=OUT/'source_comparison.json';equivalent=False
    if comparison.is_file():
     data=json.loads(comparison.read_text())
     equivalent=data.get('status')=='PASS' and any(
      v['name']==p.name and v.get('role')=='EXACT_EQUALITY'
      and v.get('historical_sha256')==item['sha256']
      and v['sha256']==sha(preferred) for v in data['files'])
    if not equivalent:return False
   elif preferred.is_file():p=preferred
  if not p.is_file() or sha(p)!=item['sha256']:return False
 return True
