"""Tie every imported generic W/h coefficient to independently rebuilt sources."""
from sage.all import *
from pathlib import Path
import json,hashlib
h=Path(__file__).absolute().parent;research=h.parents[2]
fresh=research/'full_recompute_20261009'/'rank3'/'source_replays';orig=research/'round10'
checks=[];records=[]
def use(job,name):
 root=fresh/job
 status=json.loads((root/'replay_record.json').read_text())
 assert status['status']=='PASS' and not (root/'execution.lock').exists(),job
 f=root/'generated'/'round10'/(name+'.sobj')
 records.append(dict(job=job,file=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 return load(str(f))
def same(label,a,b):
 assert a==b,label
 checks.append(label)
fw=use('pure_x_root_generic_w_reader','pure_x_root_generic_w_reader')
ow=load(str(orig/'pure_x_root_generic_w_reader.partial.sobj'))
same('all three full positive-kappa G sources used by horizontal producer',fw['full_G_sources'],ow['full_G_sources'])
ob=load(str(orig/'pure_x_root_generic_w_bridge.sobj'))
for key in ['combination','incidence_cofactor']:
 same('full W reader bridge '+key,fw[key],ob[key])
# Parent source-link report verifies the actual derivative/negative-layer
# coefficient arrays in the bridge; it is required, not a sampled witness.
pr=research/'full_recompute_20261009'/'rank3'/'fresh_remaining_source_links.json'
pd=json.loads(pr.read_text());assert pd['status']=='PASS'
records.append(dict(file=str(pr),sha256=hashlib.sha256(pr.read_bytes()).hexdigest()))
fh=use('pure_x_super_reader_generic_h_kappa_bridge','pure_x_super_reader_generic_h_kappa_bridge')
oh=load(str(orig/'pure_x_super_reader_generic_h_compact_certificate.sobj'))
for key in ['support','generic_combination','incidence_cofactor']:
 same('whole h source bridge '+key,fh[key],oh[key])
report=dict(status='PASS',checks=checks,records=records)
(h/'pure_x_external_bridges.json').write_text(json.dumps(report,indent=2))
print('PURE AXIS FULL SOURCE BRIDGES PASS',len(checks),flush=True)
