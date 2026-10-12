"""Rebuild the cyclic wall's entire 28x8 source matrix before its kernel use."""
from pathlib import Path
from sage.all import *
import json,hashlib
h=Path(__file__).absolute().parent;research=h.parents[2];o=research/'round10'
job=research/'full_recompute_20261009'/'rank3'/'source_replays'/'independent_subcritical_single_RR_conic_relative_reader'
assert json.loads((job/'replay_record.json').read_text())['status']=='PASS'
nf=job/'generated'/'round10'/'independent_subcritical_single_RR_conic_relative_reader.sobj'
d=load(str(nf));wall=load(str(o/'subcritical_single_cyclic_full_conic_wall_probe.sobj'));cyc=load(str(o/'subcritical_single_cyclic_conic_matrix_probe.sobj'));nine=load(str(o/'subcritical_single_cyclic_conic_nine_probe.sobj'))
E=wall['ring'];W=E.base_ring();K=W.base_ring();k,r=K.gens();w=E(W.gen());x=E.gen()
old=d['ring'];OK=old.base_ring();ph=OK.hom([E(k),E(r),E.one(),-E.one(),w],E)
def cv(v):return ph(v[0])+ph(v[1])*x
M=d['matrix'].apply_map(cv)
checks=[]
def ck(label,v):assert v,label;checks.append(label)
for i in range(28):
 for j in range(8):ck('full wall source matrix '+str((i,j)),M[i,j]==wall['items'][0]['matrix'][i,j])
for pos,ref in enumerate([cyc,nine]):
 CE=ref['ring'];CK=CE.base_ring();mp=CK.hom([E(k),E(r),w],E)
 def cp(v):return mp(v[0])+mp(v[1])*x
 for i in range(28):
  for j in range(8 if pos==0 else 7):ck('cyclic source range '+str((pos,i,j)),cp(ref['matrix'][i,j])==M[i,j])
B=vector(E,[0,0,0,0,0,0,-2,1]);T=wall['items'][0]['kernel'][0]
ck('actual reconstructed source B kernel',not bool(M*B))
ck('actual reconstructed source T kernel',not bool(M*T))
g=load(str(h/'replay'/'round10'/'subcritical_single_cyclic_conic_G_normal_cover.sobj'))
ck('literal G direction B',g['directions'][0]==B)
ck('literal G direction TplusB',g['directions'][1]==T+B)
# Preserve arbitrary w before either material-wall reduction.
CE=cyc['ring'];CK=CE.base_ring();ckk,crr,cww=CK.gens();xx=CE.gen()
mp=OK.hom([CE(ckk),CE(crr),CE.one(),-CE.one(),CE(cww)],CE)
for i in range(28):
 for j in range(8):ck('whole offwall cyclic matrix '+str((i,j)),mp(d['matrix'][i,j][0])+mp(d['matrix'][i,j][1])*xx==cyc['matrix'][i,j])
NE=nine['ring'];NK=NE.base_ring();nkk,nrr,nww=NK.gens();nx=NE.gen();cn=CK.hom([NE(nkk),NE(nrr),NE(nww)],NE)
for i in range(28):
 for j in range(7):ck('whole cyclic nine firstseven '+str((i,j)),cn(cyc['matrix'][i,j][0])+cn(cyc['matrix'][i,j][1])*nx==nine['matrix'][i,j])
wall2=load(str(o/'independent_subcritical_single_RR_conic_wall.sobj'));E2=wall2['ring'];W2=E2.base_ring();K2=W2.base_ring();k2,r2=K2.gens();w2=E2(W2.gen());x2=E2.gen();p2=OK.hom([E2(k2),E2(r2),E2(2),-E2.one(),w2],E2)
for i in range(28):
 for j in range(8):ck('whole H2 wall matrix '+str((i,j)),p2(d['matrix'][i,j][0])+p2(d['matrix'][i,j][1])*x2==wall2['matrix'][i,j])
for v in wall2['kernel']:ck('H2 actual normal kernel '+str(len(checks)),not bool(wall2['matrix']*vector(E2,v)))
report=dict(status='PASS',checks=checks,source=str(nf),source_sha256=hashlib.sha256(nf.read_bytes()).hexdigest(),scope='All 28x8 source entries, cyclic range and both true kernel directions. No stored determinant or rank is used as a source assertion.')
(h/'conic_wall_source_bridge.json').write_text(json.dumps(report,indent=2))
print('CONIC WALL FULL SOURCE BRIDGE PASS',len(checks),flush=True)
